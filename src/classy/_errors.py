"""The exceptions classy raises.

They are defined once, here, and both the compiled core and the package raise them, so that
`except classy.CosmoSevereError` catches whatever is raised. The C++ code throws
std::invalid_argument (class_test_severe), which becomes CosmoSevereError, and
std::runtime_error (class_test), which becomes CosmoComputationError.
"""


class CosmoError(Exception):
    def __init__(self, message=None):
        if message is None:
            message = ''
        elif isinstance(message, bytes):
            try:
                # If message has garbage inside, decode may fail.
                message = message.decode()
            except UnicodeDecodeError:
                pass

        self.message = message

    def __str__(self):
        return '\n\nError in Class: ' + self.message


class CosmoSevereError(CosmoError):
    """
    Raised when Class failed to understand one or more input parameters.

    This case would not raise any problem in Class default behaviour. However,
    for parameter extraction, one has to be sure that all input parameters were
    understood, otherwise the wrong cosmological model would be selected.
    """
    pass


class CosmoComputationError(CosmoError):
    """
    Raised when Class could not compute the cosmology at this point.

    This will be caught by the parameter extraction code to give an extremely
    unlikely value to this point
    """

    pass
