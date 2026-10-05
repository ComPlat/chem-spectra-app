import json
import logging
import os

import numpy as np
from scipy import signal

from chem_spectra.lib.converter.datatable import DatatableModel
from chem_spectra.lib.shared.calc import (to_float, cal_cyclic_volta_shift_prev_offset_at_index)
from chem_spectra.lib.converter.jcamp.data_parse import make_ni_data_ys, make_ni_data_xs
from chem_spectra.lib.converter.jcamp.techniques import technique_for
from chem_spectra.lib.converter.jcamp.records import UNIT_RECORDS
# defined in share, re-exported here: the app's error handler
# (chem_spectra/__init__.py), transformer.py and the tests all import it
# from this module, which is where most of the raises are
from chem_spectra.lib.converter.share import UnconvertibleSpectrum  # noqa: F401


# Real absorbance runs roughly 0-3; beyond this T = 10**(-A) underflows and
# the spectrum becomes a flat line.
ABSORBANCE_CEILING = 10.0

# y units that already mean transmittance, after _normalise_unit. Instrument
# exports spell it many ways -- '%T', 'T', 'TRANSMISSION', 'Transmission (%)',
# '% TRANSMITTANCE' -- and a substring test for 'TRANSMITTANCE' missed most of
# them, leaving the median heuristic, which is unreliable on heavily absorbing
# samples, as the only guard. An exact set rather than a looser substring
# match, so an unrelated unit cannot block a conversion by accident.
TRANSMITTANCE_UNITS = frozenset({
    'T', 'TRANSMITTANCE', 'TRANSMISSION',
    'PERCENTT', 'PERCENTTRANSMITTANCE', 'PERCENTTRANSMISSION',
    'TRANSMITTANCEPERCENT', 'TRANSMISSIONPERCENT',
})


def absorbance_to_percent_transmittance(values):
    """%T = 100 * 10**(-A)."""
    return 100.0 * np.power(10.0, -np.asarray(values, dtype=float))


def _normalise_unit(unit):
    """Upper case, letters and digits only: '%T' -> 'T',
    'Transmission (%)' -> 'TRANSMISSION'."""
    return ''.join(ch for ch in str(unit).upper() if ch.isalnum())


def is_transmittance_unit(unit):
    return bool(unit) and _normalise_unit(unit) in TRANSMITTANCE_UNITS


# Deliberately short. 'AU' is arbitrary units, not absorbance units, and a
# bare 'A' is as likely to be amperes; guessing either way is worse than
# falling back to the technique.
ABSORBANCE_UNITS = frozenset({'ABSORBANCE', 'ABS'})


def is_absorbance_unit(unit):
    return bool(unit) and _normalise_unit(unit) in ABSORBANCE_UNITS


logger = logging.getLogger(__name__)

data_type_json = os.path.join(os.path.dirname(__file__), 'data_type.json')

class JcampTechniqueConverter:
    """Reads a JCAMP file for any technique dispatched by its descriptor.

    Formerly JcampTechniqueConverter, where NI meant "NMR & IR" -- the two
    techniques it handled when it was written. It now covers every technique
    in SPECTRUM_TECHNIQUES except mass spectrometry, which still has its own
    converter and composer. The name states where the pipeline is going:
    MS is a technique too, and folding it in is a deferred action in
    IMPLEMENTATION-PLAN.spectrum-kind-refactor.md.
    """

    def __init__(self, base):
        self.base = base
        self.params = base.params
        self.datatypes = base.datatypes
        self.datatype = base.datatype
        self.dataclass = base.dataclass
        self.data_format = base.data_format
        self.typ = base.typ
        # resolved once, by __target_block_records: the records cannot
        # change, and the mismatch warning should be said once if at all
        self.__target_records = None
        self.target_idx = self.__index_target()
        self.dic = base.dic
        self.data = make_ni_data_ys(base, self.target_idx)
        self.title = base.title
        # the descriptor travels with the flags it backs; without it the
        # composer falls back to UNKNOWN_TECHNIQUE and draws every technique
        # with the generic-curve defaults
        self.technique = getattr(base, 'technique', None) or technique_for(self.typ)
        self.ncl = base.ncl
        self.solv_peaks = base.solv_peaks
        # - - - - - - - - - - -
        self.fname = base.fname
        # set by __to_transmittance *this run*, read by __set_label: only a
        # conversion we performed may rename the unit.
        self.converted_to_transmittance = False
        # ...whereas the record says the data *is* transmittance, whoever
        # converted it and whenever. A file we composed earlier already
        # carries it, and recomposing must not throw that away: see below.
        self.transmittance_recorded = self.__declared_flag('$CSTRANSMITTANCE')
        # A drawing instruction, carried to the composer and written into the
        # composed file as ##$CSINVERTY. It never reaches self.ys: inversion
        # is a viewport concern, and `max(y) - y` on the stored array destroys
        # the baseline, skews the area-under-curve integration and travels
        # into every downstream consumer of the exported JCAMP.
        #
        # Display only, so nothing computed from the data follows it: peak
        # picking, integration and the peak tables all run on self.ys, and
        # whether peaks are maxima or dips is decided by what the y axis
        # measures (see __peaks_point_down), never by this flag. Under #298,
        # which mirrored the array, the picker ran on the mirrored data and
        # picked the other polarity; it no longer does.
        #
        # It is set by the request *or* by the file's own record. Without the
        # second half the flag is write-only: every pass through this app
        # recomposes, and a recompose carries no `invert_y`, so the record
        # written on one request is dropped on the next. A viewing preference
        # that does not survive a round trip is not a preference. The file's
        # declaration is authoritative unless the request overrides it, which
        # is the rule #312 applied to the stored point order. Overriding works
        # both ways: `invert_y` is None when not sent, and an explicit False
        # clears the record.
        requested = base.params.get('invert_y')
        self.draw_y_inverted = (self.__declared_flag('$CSINVERTY')
                                if requested is None else bool(requested))
        self.block_count = self.__count_block()
        self.threshold = self.technique.threshold
        self.obs_freq = self.__set_obs_freq()
        self.x_unit = self.__set_x_unit()
        self.ys = self.__read_ys()
        # after __read_ys, which is where a conversion happens
        self.peaks_point_down = self.__peaks_point_down()
        self.xs = self.__read_xs(base)
        self.__check_cylic_volta_shifted_info()

        self.factor = self.__set_factor(base)
        self.__set_first_last_xs()
        self.clear = self.__refresh_solvent()
        self.boundary = self.__find_boundary()
        self.label = self.__set_label()
        self.simu_peaks = base.simu_peaks
        self.auto_peaks = None
        self.edit_peaks = None
        self.itg_table = []
        self.mpy_itg_table = []
        self.mpy_pks_table = []
        self.max_min_peaks_table = []
        self.datatable = self.__set_datatable()
        self.__read_peak_from_file()
        self.__read_integration_from_file()
        self.__read_multiplicity_from_file()
        self.__read_voltammetry_data_from_file()

    def __read_user_data_type_mapping(self):
        user_dt_mapping = self.params.get('user_data_type_mapping')
        if user_dt_mapping == '' or user_dt_mapping is None:
            return ''
        else:
            return json.loads(user_dt_mapping)['datatypes']

    def __index_target(self):
        """Index of the block holding the primary measurement.

        Only PRIMARY datatypes belong in data_type.json. Auxiliary blocks that
        sit alongside a primary one in the same file -- NMR FID, NMR PEAK
        TABLE, NMP PEAK ASSIGNMENTS, INFRARED PEAK TABLE, INFRARED
        INTERFEROGRAM -- are deliberately absent from it, because this picks
        the first RECOGNISED block and would otherwise read the auxiliary one
        instead of the spectrum. test_auxiliary_blocks_stay_unmapped enforces
        this. A datatype missing from the map is not an error: the file takes
        the generic curve path and JcampBaseConverter logs it.
        """
        if self.params.get('user_data_type_mapping'):
            data_type_mappings = self.__read_user_data_type_mapping()
            target = data_type_mappings.values()
            target_topics = [value.upper() for values in target for value in values]
        else:
            with open(data_type_json, 'r') as mapping_file:
                target = json.load(mapping_file).get("datatypes").values()
                target_topics = [value.upper() for values in target for value in values]

        # Take the first recognised block in the file's own order. The old
        # loop had no break, so the LAST entry of the flattened mapping won
        # instead -- an order nobody chose, and one that could disagree with
        # the classification in JcampBaseConverter.__set_datatype.
        idx = None
        for pos, dt in enumerate(self.datatypes):
            if dt in target_topics:
                idx = pos
                break

        if idx is None:
            # Nothing in this file is a recognised datatype. Fall back to the
            # first block and skip the LINK offset below: it would drive the
            # index negative and silently read the last block instead. The
            # unrecognised datatype is logged by JcampBaseConverter.
            #
            # `target_idx = 0` means the first *data* block, so datatype_pos
            # has to name the same one. Taking position 0 instead named the
            # outer LINK wrapper, which declares no units -- so an unmapped
            # datatype in a LINK file lost its ##XUNITS=/##YUNITS= and was
            # relabelled PPM/ARBITRARY. Every composed file is LINK-wrapped,
            # so a single-block file survived its first compose and lost its
            # units on the next one.
            self.datatype_pos = next(
                (pos for pos, dt in enumerate(self.datatypes)
                 if dt != 'LINK'), 0)
            return 0

        # The position in the file's own ##DATA TYPE= sequence, which the
        # LINK subtraction below throws away. __target_block_records needs
        # it: nmrglue's DATATYPE list is exactly the blocks that declare one,
        # in file order, so this indexes them directly.
        self.datatype_pos = idx

        if 'LINK' in self.datatypes:
            count_link = self.datatypes.count('LINK')
            idx -= count_link

        return max(idx, 0)

    def __count_block(self):
        count = 1
        try:
            count = int(self.dic['BLOCKS'][0])
        except:  # noqa
            pass

        return count

    def __read_xs(self, base):  # TBD
        if (base.data_format == '(XY..XY)'):
            xs = make_ni_data_xs(base)
            return xs

        beg_pt = None
        end_pt = None
        idx = self.target_idx

        if beg_pt is None:
            try:
                obs_freq = self.obs_freq
                shift = float(self.dic['$OFFSET'][idx])
                beg_pt = float(
                    self.dic['FIRST'][idx].replace(' ', '').split(',')[0]
                ) / obs_freq
                end_pt = float(
                    self.dic['LAST'][idx].replace(' ', '').split(',')[0]
                ) / obs_freq
                shift = beg_pt - shift
                beg_pt = beg_pt - shift
                end_pt = end_pt - shift
            except:  # noqa
                pass

        if beg_pt is None:  # MNova
            try:
                obs_freq = self.obs_freq
                beg_pt = float(
                    self.dic['FIRST'][idx].replace(' ', '').split(',')[0]
                ) / obs_freq
                end_pt = float(
                    self.dic['LAST'][idx].replace(' ', '').split(',')[0]
                ) / obs_freq
            except:  # noqa
                pass

        if beg_pt is None:
            try:
                beg_pt = to_float(self.dic['FIRSTX'][idx])
                end_pt = to_float(self.dic['LASTX'][idx])
            except:  # noqa
                pass
            
        if beg_pt is None:
            try:
                while len(self.dic['FIRSTX']) <= idx:
                    self.dic['FIRSTX'].insert(0, '')
                while len(self.dic['LASTX']) <= idx:
                    self.dic['LASTX'].insert(0, '')
                beg_pt = to_float(self.dic['FIRSTX'][idx])
                end_pt = to_float(self.dic['LASTX'][idx])
            except:  # noqa
                pass

        # Store the points the way the technique is conventionally drawn.
        # Which way that is comes from `x_reversed`; whether it may be
        # normalised at all comes from `store_in_drawn_order`, because for
        # some techniques the traversal order is data (see the registry).
        #
        # This used to ask `em_wave` for both, which is a header-shape
        # grouping and says nothing about an axis. Written for INFRARED
        # alone in 2019 (40da642), it grew to the grouping when Raman joined
        # it (af68fdd), and swept UV-Vis in when UV/VIS was first recognised
        # (387083a -- a commit about identifiers, thresholds and bin counts).
        # UV-Vis ascends by convention, so the grouping stored it backwards.
        if (self.technique.store_in_drawn_order
                and beg_pt is not None and end_pt is not None
                and beg_pt != end_pt
                and (beg_pt > end_pt) != self.technique.x_reversed):
            beg_pt, end_pt = end_pt, beg_pt
            self.ys = self.ys[::-1]

        num_pt = self.ys.shape[0]
        
        x = np.linspace(
            beg_pt + self.params['delta'],
            end_pt + self.params['delta'],
            num=num_pt,
            endpoint=True
        )

        if self.x_unit == 'HZ':
            x = x / self.obs_freq
        
        return x

    def __declared_flag(self, record):
        """Whether the file carries `##<record>=true`.

        Only the real record counts. The original-metadata dump re-emits it as
        `###$CSINVERTY= true`, which parses under a different key, so reading
        it here cannot resurrect a flag from a file that never set one.
        """
        values = self.dic.get(record) or []
        return any(str(v).strip().lower() == 'true' for v in values)

    def __read_ys(self):
        """Apply the client's processing instructions, and nothing else.

        Nothing here is inferred. Until #296's predecessor this method guessed
        from the data's shape whether an infrared spectrum was absorbance and
        mirrored it silently, while the label was decided separately from the
        declared units -- so the two could, and did, disagree. Both decisions
        now belong to whoever supplies the file.
        """
        ys = self.data
        if ys is None:
            return ys

        if self.params.get('transmittance'):
            ys = self.__to_transmittance(ys)

        return ys

    def __declared_units(self):
        """The x and y units the *file* declares for the target block, each
        None when it declares none. The one reader of the unit records: the
        transmittance guard and __set_label used to read them separately, at
        different indices, so in a LINK or multi-block file the unit the
        guard refused on and the unit the composed file was labelled with
        could disagree.

        ##XUNITS=/##YUNITS= are read first, then the JCAMP 6 ##UNITS= triple,
        which overrides them -- but only records the target block itself
        declares. Another block's triple (an interferogram's `CM, VOLTS, ...`
        beside an absorbance spectrum) must not relabel the spectrum. When
        the file could not be split into blocks, the flattened nmrglue
        lists are indexed as before. Deliberately ignores the caller's
        `axesUnits`: that is a display preference, and the guard is about
        what the data already is.

        Runs before __set_label, which is why the guard cannot read
        self.label.
        """
        records = self.__target_block_records()
        x, y = records.get('XUNITS'), records.get('YUNITS')
        try:
            # Split first, strip after: stripping the whole record turned
            # `% TRANSMITTANCE` into `%TRANSMITTANCE`, and that is what went
            # into ##YUNITS. Only the space around the separators is noise.
            triple = records['UNITS'].strip()
            # Mnova ends the record with a comma: `HZ, ARBITRARY UNITS,
            # ARBITRARY UNITS,`. Rejecting it left the label to be taken from
            # whichever other block declared XUNITS/YUNITS.
            fields = [field.strip()
                      for field in triple.rstrip(',').split(',')]
            x, y, _ = fields
        except (KeyError, AttributeError, ValueError):
            pass
        return {'x': x, 'y': y}

    def __target_block_records(self):
        """The unit records declared in the target block.

        nmrglue's DATATYPE list is exactly the blocks that declare a
        ##DATA TYPE=, in file order, so `datatype_pos` -- the position the
        target was found at in that list -- indexes them directly.

        The earlier `target_idx + count('LINK')` indexed *every* block
        instead, which assumed each one declares a datatype and that the
        LINK blocks all come first. Neither is required: an outer block that
        declares a title and nothing else shifted the whole file by one, and
        the spectrum was labelled from the block after it.

        The two lists are compared before either is trusted. If they
        disagree the file is shaped in a way neither reader anticipated, so
        the flattened lists are indexed as before and the disagreement is
        logged rather than guessed at.
        """
        if self.__target_records is not None:
            return self.__target_records
        self.__target_records = self.__resolve_target_block_records()
        return self.__target_records

    def __resolve_target_block_records(self):
        blocks = getattr(self.base, 'block_records', None) or []
        declaring = [block for block in blocks if block.get('DATATYPE')]
        sequence = [block['DATATYPE'].upper() for block in declaring]
        position = getattr(self, 'datatype_pos', None)
        if blocks and sequence == self.datatypes:
            if position is not None and 0 <= position < len(declaring):
                return declaring[position]
            return {}
        if blocks:
            logger.warning(
                'the ##DATA TYPE= sequence read from the file, %s, is not the '
                'one nmrglue reports, %s; falling back to the flattened '
                'records for %r',
                sequence, self.datatypes, self.params.get('fname'),
            )

        def per_block(records):
            for idx in (self.target_idx, 0):
                try:
                    return records[idx]
                except (IndexError, KeyError, TypeError):
                    continue
            return None

        found = {key: per_block(self.dic.get(key)) for key in UNIT_RECORDS}
        return {key: value for key, value in found.items() if value}

    def __peaks_point_down(self):
        """Whether the bands of interest are dips in the stored trace.

        The quantity decides it, not the technique. `peaks_inverted` is set
        per technique because an infrared spectrum is nearly always %T, but
        that is a habit of the format, not a property of infrared: an IR
        file in absorbance has maxima, and a UV/VIS spectrum converted to %T
        has dips. Reading the technique alone meant a converted UV/VIS
        spectrum came back with its auto peaks on the baseline *between* the
        bands, and an IR absorbance file had no request that would find its
        bands at all.

        The technique stays as the fallback, for the files -- most NMR among
        them -- whose y unit says nothing about direction.
        """
        if self.converted_to_transmittance or self.transmittance_recorded:
            return True
        declared = self.__declared_units()['y']
        if is_transmittance_unit(declared):
            return True
        if is_absorbance_unit(declared):
            return False
        return self.technique.peaks_inverted

    def __peak_threshold(self):
        """`technique.threshold` is a fraction of the maximum, read against
        the technique's own polarity: 0.93 for infrared means "dips below
        93% of the maximum", 0.05 for UV/VIS means "maxima above 5% of it".
        Where the polarity we need is the other one, so is the fraction --
        otherwise a converted UV/VIS trace is searched for dips below 5% of
        the maximum, which is the floor, and nothing is found.
        """
        if self.peaks_point_down == self.technique.peaks_inverted:
            return self.threshold
        return 1.0 - self.threshold

    def __to_transmittance(self, ys):
        """T = 10**(-A). Refuses rather than returning a ruined spectrum.

        The conversion is only meaningful for real absorbance, which runs
        roughly 0-3. Asked to convert anything else it would silently produce
        a flat line, or infinities, labelled TRANSMITTANCE -- worse than any
        of the defects this change fixes.

        Every refusal is checked on the *input*, so the reason names what is
        wrong with the data rather than what went wrong arithmetically. The
        finiteness check on the result is a backstop for anything the input
        checks do not anticipate.
        """
        # Our own record is checked first: a file we converted earlier is
        # transmittance whatever unit a later recompose wrote over it.
        if self.transmittance_recorded:
            raise UnconvertibleSpectrum(
                'the file records that it was already converted to '
                'transmittance; there is nothing to convert'
            )
        declared = self.__declared_units()['y']
        if is_transmittance_unit(declared):
            # What the file says outranks what its shape suggests: this is a
            # fact, where the median test below is an inference. It is also
            # checked before the integrals: there is no point telling someone
            # to remove integrals from a file that is going to be refused as
            # already transmittance whatever they do.
            raise UnconvertibleSpectrum(
                'the file already declares its y axis as {!r}; there is '
                'nothing to convert'.format(declared)
            )
        if self.__carries_integrals():
            # An area under absorbance is proportional to concentration; the
            # same region of the %T trace has no such meaning, and %T is not
            # linear in A, so the areas cannot be carried across either.
            raise UnconvertibleSpectrum(
                'the spectrum carries integrals or multiplets, which have no '
                'meaning in transmittance; remove them before converting'
            )

        ys = self.__refuse_unless_absorbance(ys)
        y_max = float(np.max(ys))
        y_min = float(np.min(ys))

        if float(np.median(ys)) >= 0.5 * y_max:
            raise UnconvertibleSpectrum(
                'already appears to be transmittance (baseline near the '
                'maximum); there is nothing to convert'
            )
        # Percent, not the 0-1 ratio: %T is how instruments commonly present
        # transmittance, and absorbance itself runs 0-2.5, so a 0-1 array is
        # routinely misread as absorbance. '% TRANSMITTANCE' is NOT a JCAMP-DX
        # unit, though: 4.24 (6.2.2) lists TRANSMITTANCE only as the ratio
        # I_T/I_0, beside REFLECTANCE, ABSORBANCE, KUBELKA-MUNK and ARBITRARY
        # UNITS. A strict reader will not recognise it as transmittance, and
        # one that maps it to TRANSMITTANCE sees values 100x too large.
        transmittance = self.__to_percent(ys)

        self.converted_to_transmittance = True
        self.transmittance_recorded = True
        return transmittance

    @staticmethod
    def __refuse_unless_absorbance(values, what=None):
        """The range checks both the trace and the peak table need.

        Every refusal names what is wrong with the *input*, so the reason
        says what the data is rather than what the arithmetic did. `what`
        names the table when it is not the trace; the trace's own wording is
        unchanged, because it is what the API has been answering with.
        """
        where = '{}: '.format(what) if what else ''
        values = np.asarray(values, dtype=float)
        if values.size == 0:
            return values
        if not np.isfinite(values).all():
            raise UnconvertibleSpectrum(
                where + 'the series contains non-finite values, so it cannot '
                'be absorbance'
            )
        y_max = float(np.max(values))
        if y_max > ABSORBANCE_CEILING:
            raise UnconvertibleSpectrum(
                where + 'y values up to {:g} are not absorbance; '
                'transmittance would underflow to zero'.format(y_max)
            )
        # Absorbance dips slightly below zero from baseline drift, but not
        # far: A = -10 already means T = 10**10, which is not a
        # transmittance.
        y_min = float(np.min(values))
        if y_min < -ABSORBANCE_CEILING:
            raise UnconvertibleSpectrum(
                where + 'y values down to {:g} are not absorbance; '
                'transmittance would overflow'.format(y_min)
            )
        return values

    @classmethod
    def __to_percent(cls, values, what=None):
        """Guarded conversion. The finiteness check on the result is a
        backstop for anything the input checks do not anticipate."""
        where = '{}: '.format(what) if what else ''
        values = cls.__refuse_unless_absorbance(values, what)
        converted = absorbance_to_percent_transmittance(values)
        if converted.size and not np.isfinite(converted).all():
            raise UnconvertibleSpectrum(
                where + 'the conversion produced non-finite values; the '
                'series is not absorbance'
            )
        return converted

    def __carries_integrals(self):
        """Integrals or multiplets the composed file would carry.

        Not "present anywhere": a request that clears the table clears it.
        The composer writes nothing when an edited table arrives empty --
        `gen_integration_info` for the integrals, `gen_mpy_integ_info` and
        `gen_mpy_peaks_info` for the multiplets -- so refusing on the file's
        stale record would refuse a conversion over a table on its way out.

        `edited` absent is not `edited` false: parse_params supplies a
        default with no such key, so a request that simply does not mention
        integrals leaves the file's record standing, as it should.
        """
        for param, record in (('integration', '$OBSERVEDINTEGRALS'),
                              ('multiplicity', '$OBSERVEDMULTIPLETS')):
            sent = self.params.get(param) or {}
            if sent.get('stack'):
                return True
            if sent.get('edited'):
                # an explicit, empty edit: the table is being removed
                continue
            if self.__record_has_rows(record):
                return True
        return False

    def __record_has_rows(self, record):
        """Whether a peak-table record holds any data rows.

        By shape, not by position. `$OBSERVEDINTEGRALS` opens with an
        `(X Y Z)` header and `$OBSERVEDMULTIPLETS` has none at all, so
        skipping the first line read a one-row multiplet table as empty.
        A data row is parenthesised and carries at least one digit, which
        no column header does.
        """
        for value in self.dic.get(record) or []:
            for line in str(value).split('\n'):
                line = line.strip()
                if (line.startswith('(') and
                        any(char.isdigit() for char in line)):
                    return True
        return False

    def __find_boundary(self):
        return {
            'x': {
                'max': self.xs.max(),
                'min': self.xs.min(),
            },
            'y': {
                'max': self.ys.max(),
                'min': self.ys.min(),
            },
        }

    def __set_label(self):
        declared = self.__declared_units()
        x, y = declared['x'] or 'PPM', declared['y'] or 'ARBITRARY'
        target = {
            'x': 'PPM' if x.upper() == 'HZ' else x,
            # Bruker LINK files already arrive space-stripped, so this
            # compared the squeezed spelling. Now that the triple keeps its
            # spaces, the comparison has to do the squeezing itself or
            # Mnova's `ARBITRARY UNITS` stops matching and the two sources
            # label the same quantity differently again.
            'y': ('ARBITRARY' if y.upper().replace(' ', '') == 'ARBITRARYUNITS'
                  else y),
        }

        if self.technique.x_axis == 'xrd':
            target['x'] = '2Theta'
            
        if 'axesUnits' in self.params and self.params['axesUnits'] is not None:
          axesUnits = self.params['axesUnits']
          xUnit, yUnit = axesUnits['xUnit'], axesUnits['yUnit']
          if xUnit != '':
            target['x'] = xUnit
          if yUnit != '':
            target['y'] = yUnit

        # A conversion is a fact, so it outranks axesUnits, which is a
        # preference -- whether it happened on this request or on an earlier
        # one that left ##$CSTRANSMITTANCE behind. Checking only this run let
        # a recompose label %T data with the caller's absorbance unit while
        # still writing the record.
        if self.converted_to_transmittance or self.transmittance_recorded:
            target['y'] = '% TRANSMITTANCE'

        return target

    def __set_obs_freq(self):
        obs_freq = None
        try:
            obs_freq = float(self.dic['.OBSERVEFREQUENCY'][self.target_idx])
        except:  # noqa
            try:
                 obs_freq = float(self.dic['.OBSERVEFREQUENCY'][0])
            except:  # noqa
              pass
        try:
            if obs_freq is None:
                obs_freq = float(self.dic['$SFO1'][0])
        except:  # noqa
            pass

        return obs_freq

    def __set_factor(self, base):
        factor = {
            'x': 1.0,
            'y': 1.0,
        }

        if (self.data_format and self.data_format == '(XY..XY)'):
            return factor

        try:
            factor = {
                'x': to_float(self.dic['XFACTOR'][0]),
                'y': to_float(self.dic['YFACTOR'][0]),
            }
        except:  # noqa
            try:
                factor_line = self.dic['FACTOR']
                real_factor = factor_line[0].split(",")
                factor = {
                    'x': to_float(real_factor[0]),
                    'y': to_float(real_factor[1]),
                }
            except:
                pass
        
        if factor['y'] == 1.0 and not isinstance(self.base.data, dict):
            factor['y'] = self.data.max() / 1000000.0

        return factor

    def __set_x_unit(self):
        x_unit = None
        
        if 'axesUnits' in self.params and self.params['axesUnits'] is not None:
          axesUnits = self.params['axesUnits']
          xUnit = axesUnits['xUnit']
          if xUnit != '':
            return xUnit

        try: # jcamp version 6
            units = self.dic['UNITS']
            array_unit = units[0].split(',')
            x_unit = (array_unit[0].upper()).strip()
        except: # noqa
            pass

        if (x_unit is None):
            try:
                x_unit = self.dic['XUNITS'][self.target_idx].upper()
            except:  # noqa
                try:
                     x_unit = self.dic['XUNITS'][0].upper()
                except:
                    pass

        return x_unit

    def __read_auto_peaks(self):
        if self.params['clear'] or self.clear:
            return

        try:  # legacy
            auto_x = []
            auto_y = []
            pas = self.dic['PEAKASSIGNMENTS'][0].split('\n')[1:]
            for pa in pas:
                info = pa.replace('(', '').replace(')', '') \
                            .replace(' ', '').split(',')
                auto_x.append(float(info[0]))
                auto_y.append(float(info[1]))
            if len(auto_x) == 0:
                return
            self.auto_peaks = {'x': auto_x, 'y': auto_y}
        except:  # noqa
            pass

        try:  # mnova
            if self.auto_peaks is None:
                auto_x = []
                auto_y = []
                if len(self.dic['PEAKTABLE']) == 0:
                    return
                pas = self.dic['PEAKTABLE'][1].split('\n')[1:]
                for pa in pas:
                    info = pa.replace(' ', '').split(',')
                    auto_x.append(float(info[0]))
                    auto_y.append(float(info[1]))
                if len(auto_x) == 0:
                    return
                self.auto_peaks = {'x': auto_x, 'y': auto_y}
        except:  # noqa
            pass

    def __read_edit_peaks(self):
        if self.params['clear'] or self.clear:
            return

        try:  # legacy
            edit_x = []
            edit_y = []
            pas = self.dic['PEAKASSIGNMENTS'][1].split('\n')[1:]
            for pa in pas:
                info = pa.replace('(', '').replace(')', '') \
                            .replace(' ', '').split(',')
                edit_x.append(float(info[0]))
                edit_y.append(float(info[1]))
            if len(edit_x) == 0:
                return
            self.edit_peaks = {'x': edit_x, 'y': edit_y}
        except:  # noqa
            pass

        try:  # mnova
            if self.edit_peaks is None:
                edit_x = []
                edit_y = []
                pas = self.dic['PEAKTABLE'][0].split('\n')[1:]
                for pa in pas:
                    info = pa.replace(' ', '').split(',')
                    edit_x.append(float(info[0]))
                    edit_y.append(float(info[1]))
                if len(edit_x) == 0:
                    return
                self.edit_peaks = {'x': edit_x, 'y': edit_y}
        except:  # noqa
            pass

    def __parse_edit(self):
        peaks_str = self.params['peaks_str']
        if not peaks_str:
            self.edit_peaks = {'x': [], 'y': []}
            return
        edit_x = []
        edit_y = []
        for p in peaks_str.split('#'):
            info = p.split(',')
            edit_x.append(float(info[0]))
            edit_y.append(float(info[1]))
        self.edit_peaks = {'x': edit_x, 'y': edit_y}

    def __exec_peak_picking_logic(self, refresh_solvent=False):
        # Polarity comes from what the y axis measures (__peaks_point_down)
        # and, for the second pass, from `negative_peaks`. Never from
        # draw_y_inverted: invert_y is a viewport flip and the picker sees
        # the data as stored.
        max_y = np.max(self.ys)
        height = (0.2 * max_y if refresh_solvent
                  else self.__peak_threshold() * max_y)

        corr_data_ys = self.ys
        corr_height = height
        if self.peaks_point_down:
            corr_data_ys = 1 - self.ys
            corr_height = 1 - height

        peak_idxs = signal.find_peaks(corr_data_ys, height=corr_height)[0]

        min_y = np.min(self.ys)
        if self.technique.negative_peaks and (max_y * 0.4 < -min_y):
            dept_corr_data_ys = 1 - self.ys
            dept_corr_height = height
            dept_peak_idxs = signal.find_peaks(dept_corr_data_ys, height=dept_corr_height)[0]
            peak_idxs = np.unique(np.concatenate((peak_idxs, dept_peak_idxs)))
        return peak_idxs

    def __run_auto_pick_peak(self):
        peak_idxs = self.__exec_peak_picking_logic()
        auto_peaks = [{'x': self.xs[idx], 'y': self.ys[idx]} for idx in peak_idxs]
        auto_peaks.sort(key=lambda d: d['y'], reverse=True)

        if self.peaks_point_down:
            # sorted by descending y, so the deepest dips are at the end
            auto_peaks = auto_peaks[-100:]
        elif self.ncl == '13C':
            simu_length = len(self.simu_peaks)
            simu_length = simu_length if simu_length > 1 else 50
            auto_peaks = auto_peaks[:200]
            # rm solvent peaks
            edit_non_solv_peaks = []
            for peak in auto_peaks:
                not_solvent = True
                for u, v in self.solv_peaks:
                    if u < peak['x'] < v:
                        not_solvent = False
                if not_solvent:
                    edit_non_solv_peaks.append(peak)
            # as - 26.90 (range: 26.80 - 26.95)
            # bs, cs - 207.1 (207.0-207.2) + 30.9 (30.8 - 31.0)
            # ds, es, fs, gs - 60.4 (60.3 - 60.5) + 14.2 (14.1 - 14.3) + 21.1 (20.9 - 21.2) + 171.3 (171.2-171.4)
            # rm impurity peaks
            imp_as, imp_bs, imp_cs, imp_ds, imp_es, imp_fs, imp_gs, edit_pure_peaks = [], [], [], [], [], [], [], []
            i, capacity, l = 0, simu_length + 10, len(edit_non_solv_peaks)
            while i < capacity and i < l:
                target = edit_non_solv_peaks[i]
                if 26.80 <= target['x'] <= 27.0:
                    imp_as.append(target)
                elif 207.0 <= target['x'] <= 207.2:
                    imp_bs.append(target)
                elif 30.8 <= target['x'] <= 31.0:
                    imp_cs.append(target)
                elif 60.3 <= target['x'] <= 60.5:
                    imp_ds.append(target)
                elif 14.1 <= target['x'] <= 14.3:
                    imp_es.append(target)
                elif 20.9 <= target['x'] <= 21.2:
                    imp_fs.append(target)
                elif 171.1 <= target['x'] <= 171.4:
                    imp_gs.append(target)
                else:
                    edit_pure_peaks.append(target)
                    i += 1
                    continue
                i += 1
                capacity += 1
            if not (imp_bs and imp_cs):
                edit_pure_peaks = edit_pure_peaks + imp_bs + imp_cs
            if not (imp_ds and imp_es and imp_fs and imp_gs):
                edit_pure_peaks = edit_pure_peaks + imp_ds + imp_es + imp_fs + imp_gs
            edit_pure_peaks.sort(key=lambda d: d['y'], reverse=True)
            edit_peaks = edit_pure_peaks[:simu_length]
            # assign to edit_peaks
            edit_x = [peak['x'] for peak in edit_peaks]
            edit_y = [peak['y'] for peak in edit_peaks]
            self.edit_peaks = {'x': edit_x, 'y': edit_y}
            auto_peaks = auto_peaks[:100]
        elif self.ncl == '1H':
            auto_peaks = auto_peaks[:100]
            edit_non_solv_peaks = []
            for peak in auto_peaks:
                not_solvent = True
                for u, v in self.solv_peaks:
                    if u < peak['x'] < v:
                        not_solvent = False
                if not_solvent:
                    edit_non_solv_peaks.append(peak)
            edit_x = [peak['x'] for peak in edit_non_solv_peaks]
            edit_y = [peak['y'] for peak in edit_non_solv_peaks]
            self.edit_peaks = {'x': edit_x, 'y': edit_y}
        else:
            auto_peaks = auto_peaks[:100]

        auto_x = [peak['x'] for peak in auto_peaks]
        auto_y = [peak['y'] for peak in auto_peaks]
        self.auto_peaks = {'x': auto_x, 'y': auto_y}

    def __set_datatable(self):
        y_factor = self.factor and self.factor['y']
        y_factor = y_factor or 1.0
        if (self.data_format and self.data_format == '(XY..XY)'):
            return DatatableModel().encode(
                self.ys,
                y_factor,
                self.xs,
                True
            )
        return DatatableModel().encode(
            self.ys,
            y_factor
        )

    def __read_peak_from_file(self):
        self.__read_auto_peaks()
        self.__read_edit_peaks()
        if self.converted_to_transmittance:
            # Every peak read so far was picked on the absorbance trace. The
            # automatic ones are re-picked, because an absorbance band and a
            # %T dip are not at the same place; the old table is wrong in
            # position as well as scale.
            self.auto_peaks = None
        if not self.auto_peaks or not self.params['delta'] == 0.0:
            self.__run_auto_pick_peak()
        if self.params['peaks_str'] is not None:
            self.__parse_edit()
        if self.converted_to_transmittance:
            # Edited peaks are the user's choice, so they keep their x and
            # only change units -- whether they came from the file or from
            # this request's peaks_str, which was sent by an editor showing
            # the absorbance trace. Converted once, and only the table that
            # survives: converting a stored table the request is about to
            # replace could refuse the whole conversion over numbers that
            # were on their way out.
            self.edit_peaks = self.__peaks_to_transmittance(self.edit_peaks)

    def __peaks_to_transmittance(self, peaks):
        """The same guards the trace gets. A stored table that is not
        absorbance -- y = -400, or 75 -- became inf or 1e-73 and was written
        into ##PEAKTABLE as a coordinate."""
        if not peaks or not peaks.get('y'):
            return peaks
        return {
            'x': peaks['x'],
            'y': self.__to_percent(peaks['y'], 'the stored peak table').tolist(),
        }

    def __read_voltammetry_data_from_file(self):
        target = self.dic.get('$CSCYCLICVOLTAMMETRYDATA')
        if target:
            target = target[0].split('\n')
            if (len(target) > 0):
                for item in target:
                    splitted_item = item.replace('(', '').replace(')', '')
                    splitted_item = [x.strip() for x in splitted_item.split(',')]
                    splitted_item = [float(x) if x != '' else x for x in splitted_item]
                    max_peak = {'x': splitted_item[0], 'y': splitted_item[1]}
                    min_peak = {'x': splitted_item[2], 'y': splitted_item[3]}
                    pecker = {'x': splitted_item[6], 'y': splitted_item[7]}
                    if pecker['x'] != '':
                        self.max_min_peaks_table.append({'max': max_peak, 'min': min_peak, 'pecker': pecker})
                    else:   
                        self.max_min_peaks_table.append({'max': max_peak, 'min': min_peak})


    def __read_integration_from_file(self):
        target = self.dic.get('$OBSERVEDINTEGRALS')
        if target:
            self.itg_table = ['\n'.join(target[0].split('\n')[1:]), '\n']

    def __read_multiplicity_from_file(self):
        target1 = self.dic.get('$OBSERVEDMULTIPLETS')
        if target1:
            self.mpy_itg_table = target1
            self.mpy_itg_table.append('\n')
        target2 = self.dic.get('$OBSERVEDMULTIPLETSPEAKS')
        if target2:
            self.mpy_pks_table = target2
            self.mpy_pks_table.append('\n')

    #     if self.ncl == '13C' and len(self.mpy_itg_table) == 0 and len(self.mpy_pks_table) == 0:
    #         self.__add_13C_mpy_programmatically()

    # def __add_13C_mpy_programmatically(self):
    #     num_edit_peaks = len(self.edit_peaks['x'])
    #     if num_edit_peaks > 50:
    #         return

    #     str_mpy_itg = ''
    #     for idx in range(num_edit_peaks):
    #         str_mpy_itg += '({}, {}, {}, {}, 1.0, {}, s, {})\n'.format(
    #             idx + 1,
    #             self.edit_peaks['x'][idx] - 1.0,
    #             self.edit_peaks['x'][idx] + 1.0,
    #             self.edit_peaks['x'][idx],
    #             idx + 1,
    #             'ABCDEFGHIJKLMNOPQRSTUVWXYZ'[idx%26]
    #         )
    #     if str_mpy_itg:
    #         self.mpy_itg_table = [str_mpy_itg]

    #     str_mpy_pks = ''
    #     for idx in range(num_edit_peaks):
    #         str_mpy_pks += '({}, {}, {})\n'.format(
    #             idx + 1,
    #             self.edit_peaks['x'][idx],
    #             self.edit_peaks['y'][idx]
    #         )
    #     if str_mpy_pks:
    #         self.mpy_pks_table = [str_mpy_pks]

    def __refresh_solvent(self):
        if self.ncl == '13C':
            ref_name = (
                self.params['ref_name'] or
                self.dic.get('$CSSOLVENTNAME', [''])[0]
            )
            if ref_name and ref_name != '- - -':
                return
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            peak_idxs = self.__exec_peak_picking_logic(refresh_solvent=True)[:100]  # noqa: E501
            auto_peaks = [{'x': self.xs[idx], 'y': self.ys[idx], 'idx': idx} for idx in peak_idxs]  # noqa: E501
            auto_peaks.sort(key=lambda d: d['y'], reverse=True)
            # is acetone
            left, right = auto_peaks[0], auto_peaks[1]
            if left['x'] < right['x']: left, right = right, left    # noqa: E701
            diff = abs(left['x'] - right['x'])
            if 175 < diff < 177:
                self.dic['$CSSOLVENTNAME'] = ['Acetone-d6 (sep)']
                self.dic['$CSSOLVENTVALUE'] = ['29.920']
                self.dic['$CSSOLVENTX'] = ['0']
                self.solv_peaks = [(27.0, 33.0), (203.7, 209.7)]
                shift = 29.920 - right['x']
                self.xs = self.xs + shift
                return True  # self.clear
            # is chloroform
            for hpk in auto_peaks[:10]:
                x_c = hpk['x']
                peaks = [p for p in auto_peaks if x_c - 2.0 < p['x'] < x_c + 2.0]   # noqa: E501
                if len(peaks) == 3:
                    pxs = sorted(map(lambda p: p['x'], peaks))
                    diff_one = abs(pxs[0] - pxs[1])
                    diff_two = abs(pxs[1] - pxs[2])
                    if 0.2 < diff_one < 0.6 and 0.2 < diff_two < 0.6:
                        self.dic['$CSSOLVENTNAME'] = ['Chloroform-d (t)']
                        self.dic['$CSSOLVENTVALUE'] = ['77.16']
                        self.dic['$CSSOLVENTX'] = ['0']
                        self.solv_peaks = [(74.0, 80.0)]
                        shift = 77.16 - pxs[1]
                        self.xs = self.xs + shift
                        return True  # self.clear

        return False  # self.clear

    def __set_first_last_xs(self):
        self.first_x = self.xs[0]
        self.last_x = self.xs[-1]

    def __check_cylic_volta_shifted_info(self):
        if not self.technique.cyclic_voltammetry:
            return
        
        cyclicvolta_data = self.params['cyclicvolta']
        current_jcamp_idx = self.params['jcamp_idx']
        offset = cal_cyclic_volta_shift_prev_offset_at_index(cyclicvolta_data, current_jcamp_idx)
        self.xs = np.array([x - offset for x in self.xs])
