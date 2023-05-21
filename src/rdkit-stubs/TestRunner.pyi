"""
rdkit.TestRunner module¶
"""
from _typeshed import Incomplete
from rdkit import RDConfig as RDConfig

class OutputRedirectC(object):
    """
    Context manager which uses low-level file descriptors to suppress
    output to stdout/stderr, optionally redirecting to the named file(s).
    Suppress all output
    with Silence():

    <code>

    Redirect stdout to file
    with OutputRedirectC(stdout=’output.txt’, mode=’w’):

    <code>

    Redirect stderr to file
    with OutputRedirectC(stderr=’output.txt’, mode=’a’):

    <code>

    http://code.activestate.com/recipes/577564-context-manager-for-low-level-redirection-of-stdou/
    >>>"""

    ...
    outfiles: Incomplete
    combine: Incomplete
    mode: Incomplete
    saved_streams: Incomplete
    fds: Incomplete
    saved_fds: Incomplete
    null_fds: Incomplete
    null_streams: Incomplete

    def __init__(self, stdout=..., stderr=..., mode: str = ...) -> None: ...
    def __enter__(self) -> None: ...
    def __exit__(self, *args): ...

class redirect_stderr(_RedirectStream):
    """
    Context manager for temporarily redirecting stderr to another file."""

    ...
    ...

class redirect_stdout(_RedirectStream):
    """
    Context manager for temporarily redirecting stdout to another file.
    # How to send help() to stderr
    with redirect_stdout(sys.stderr):

    help(dir)

    # How to write help() to a file
    with open(‘help.txt’, ‘w’) as f:

    with redirect_stdout(f):help(pow)"""

    ...
    ...

def ReportResults(self, script, failedTests, nTests, runTime, verbose, dest): ...
def RunScript(self, script, doLongTests, verbose): ...
def RunTest(self, exeName, args, extras): ...

TEST_FAILED: int
TEST_PASSED: int
BUILD_TYPE_ENVVAR: str

def isDebugBuild(self): ...
def RunTest(self, exeName, args, extras): ...
def RunScript(self, script, doLongTests, verbose): ...
def ReportResults(self, script, failedTests, nTests, runTime, verbose, dest): ...

class _RedirectStream:
    def __init__(self, new_target) -> None: ...
    def __enter__(self): ...
    def __exit__(self, exctype, excinst, exctb) -> None: ...

class redirect_stdout(_RedirectStream):
    """
    Context manager for temporarily redirecting stdout to another file.
    # How to send help() to stderr
    with redirect_stdout(sys.stderr):

    help(dir)

    # How to write help() to a file
    with open(‘help.txt’, ‘w’) as f:

    with redirect_stdout(f):help(pow)"""

    ...
    ...

class redirect_stderr(_RedirectStream):
    """
    Context manager for temporarily redirecting stderr to another file."""

    ...
    ...

class OutputRedirectC(object):
    """
    Context manager which uses low-level file descriptors to suppress
    output to stdout/stderr, optionally redirecting to the named file(s).
    Suppress all output
    with Silence():

    <code>

    Redirect stdout to file
    with OutputRedirectC(stdout=’output.txt’, mode=’w’):

    <code>

    Redirect stderr to file
    with OutputRedirectC(stderr=’output.txt’, mode=’a’):

    <code>

    http://code.activestate.com/recipes/577564-context-manager-for-low-level-redirection-of-stdou/
    >>>"""

    ...
    outfiles: Incomplete
    combine: Incomplete
    mode: Incomplete
    saved_streams: Incomplete
    fds: Incomplete
    saved_fds: Incomplete
    null_fds: Incomplete
    null_streams: Incomplete

    def __init__(self, stdout=..., stderr=..., mode: str = ...) -> None: ...
    def __enter__(self) -> None: ...
    def __exit__(self, *args): ...
