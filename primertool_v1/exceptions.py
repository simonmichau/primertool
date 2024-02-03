# AUTHOR: Daniela Dey
# Exceptions for Primertool


class PrimertoolError(Exception):
    pass


class PrimertoolGenomeError(PrimertoolError):
    pass


class PrimertoolInputError(PrimertoolError):
    pass


class PrimertoolHGVSError(PrimertoolError):
    pass


class PrimertoolMutalyzerError(PrimertoolError):
    pass
