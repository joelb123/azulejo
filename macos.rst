MacOS Installation
==================

We recommend installation using XCode Command Line Tools.  This avoids
the heavy-weight installations of either the full XCode or Anaconda Python,
while ensuring that you will get the latest updates via XCode.

In a terminal window execute this command::

        xcode-select --install

When this command completes, edit your ``.zshrc`` (or ``.bashrc`` if
you selected the non-default shell) and enter the following::

        export PATH="/Library/Developer/CommandLineTools/usr/bin:${PATH}"
        export CFLAGS="-I/Library/Developer/CommandLineTools/Library/Frameworks/Python3.framework/Versions/Current/include/"
        export LDFLAGS="-L/Library/Developer/CommandLineTools/Library/Frameworks/Python3.framework/Versions/Current/lib"

This will ensure that you get the proper version of python3 (the one
in ``/usr/bin`` is slightly impaired, and that you have the proper
headers and libraries for compiling new code.

You will also need a ``cmake`` executable if you are to build
``MMseqs```.  Visit the `cmake download site <https://cmake.org/download/>`_ and
select the latest ``dmg`` and go through the install process.
Open the ``cmake`` application and click on the security settings
that allow it to run.  Then edit your ``.zshrc`` file again and
add the following::

        export PATH="${PATH}:/Applications/Cmake.app/Contents/bin"

Restart your shell and you should be ready to proceed with the
``azulejo`` installation instructions.

