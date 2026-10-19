"""One shape for every refusal: 422 with `{"error": "..."}`.

Before this, a file the app would not process came back as Flask's HTML error
page — 400, 403, 404 or, most often, an unhandled 500. The ELN's
`Jcamp::Create.spectrum` falls back from the `x-extra-info-json` header to
`JSON.parse(rsp.body)` and raises `json_rsp['error']`; it reaches
*"Chemspectra response missing metadata header"* only when the body is not
JSON. So every refusal reached the user as that sentence, whatever the actual
reason was.

Nothing on the ELN side branches on the status code — every check is
`code == 200` — so one status is enough and the body carries the meaning.
422 is the same status `UnconvertibleSpectrum` already returns, so there is one
convention rather than two.
"""

from flask import jsonify

REFUSED = 422


def refuse(reason):
    """Refuse the request, saying why.

    `reason` is shown to the user through the ELN, so it should name what was
    wrong with *this upload* rather than where in the code it was noticed.
    """
    return jsonify(error=reason), REFUSED


def uploaded_file(request, field='file'):
    """The upload under `field`, or None when there is nothing usable.

    `request.files['file']` raises a `BadRequestKeyError` when the field is
    absent, which Flask renders as an HTML 400. And `FileContainer` defines no
    `__bool__`, so the `if file:` guards at the call sites were always true and
    their else-branches unreachable — an upload with no filename went on to
    fail somewhere further in, as a 500.

    werkzeug's `FileStorage` is falsy when its filename is empty, so this one
    check covers both the missing field and the empty upload.
    """
    uploaded = request.files.get(field)
    return uploaded or None
