class PrimertoolError(Exception):
    pass


class PrimertoolGenomeError(PrimertoolError):
    pass


class PrimertoolInputError(PrimertoolError):
    pass


class PrimertoolMutalyzerError(PrimertoolError):
    pass
