# -*- coding: utf-8 -*-
"""Install binary dependencies."""
# standard library imports
import os
import pkg_resources
import pkgutil
import shutil
import sys
import tempfile
from pathlib import Path
from packaging import version

# third-party imports
import sh
from loguru import logger
from progressbar import DataTransferBar
from requests_download import ProgressTracker
from requests_download import download as request_download


def default_version_splitter(instring):
    """Split the version string out of version output."""
    return instring.split()[-1]


def walk_package(pkgname, root):
    """Walk through a package_resource."""
    dirs = []
    files = []
    for name in pkg_resources.resource_listdir(pkgname, root):
        fullname = root + "/" + name
        if pkg_resources.resource_isdir(pkgname, fullname):
            dirs.append(fullname)
        else:
            files.append(name)
    for new_path in dirs:
        yield from walk_package(pkgname, new_path)
    yield root, dirs, files


class DependencyInstaller(object):
    """Install and check binary dependencies."""

    def __init__(
        self,
        dependency_dict,
        pkg_name=__name__,
        install_path=None,
        bin_path=None,
        force=False,
    ):
        """Initialize dictionary of dependencies."""
        self.dependency_dict = dependency_dict
        self.force = force
        self.pkg_name = pkg_name
        self.dependencies = tuple(dependency_dict.keys())
        self.versions_checked = False
        if install_path is None:
            self.install_path = Path(sys.executable).parent.parent
        else:
            self.install_path = install_path
        if bin_path is None:
            self.bin_path = self.install_path / "bin"
        else:
            self.bin_path = bin_path
        self.bin_path_exists = self.bin_path.exists()
        self.bin_path_writable = os.access(self.bin_path, os.W_OK)
        self.bin_path_in_path = str(self.bin_path) in os.environ["PATH"].split(
            ":"
        )

    def check_all(self):
        """Check all depenencies for existence and version."""
        for dep in self.dependencies:
            target_version = version.parse(
                self.dependency_dict[dep]["version"]
            )
            version_command = self.dependency_dict[dep]["version_command"]
            self.dependency_dict[dep]["installed"] = not self.force
            for bin in self.dependency_dict[dep]["binaries"]:
                if sh.which(bin) == None:
                    logger.debug(
                        f"Binary {bin} of dependency {dep} is not" " installed"
                    )
                    self.dependency_dict[dep]["installed"] = False
                else:
                    exe = sh.Command(bin)
                    ver_out = exe(*version_command, _err_to_out=True).rstrip(
                        "\n"
                    )
                    if "version_parser" in self.dependency_dict[dep]:
                        installed_version = version.parse(
                            self.dependency_dict[dep]["version_parser"](
                                ver_out
                            )
                        )
                    else:
                        installed_version = version.parse(
                            default_version_splitter(ver_out)
                        )
                    if installed_version == target_version:
                        ver_str = (
                            f"{bin} version at recommended version"
                            f" {installed_version}"
                        )
                    elif installed_version < target_version:
                        ver_str = (
                            f"{bin} installed {installed_version} <  target"
                            f" {target_version}."
                        )
                        self.dependency_dict[dep]["installed"] = False
                    elif installed_version > target_version:
                        ver_str = (
                            f"installed {installed_version} exceeds target"
                            f" {target_version}."
                        )
                    print(f"{dep}: {exe} {ver_str}")
        self.versions_checked = True
        # Check that bin directory exists and is writable.
        if self.bin_path_exists:
            bin_path_state = "exists, "
        else:
            bin_path_state = "doesn't exist, "
        if self.bin_path_writable:
            bin_path_state += "writable, "
        else:
            bin_path_state += "not writable, "
        if self.bin_path_in_path:
            bin_path_state += "in path."
        else:
            bin_path_state += "not in path."
            logger.debug(f"Bin dir '{self.bin_path}' {bin_path_state}")

    def install_list(self, deplist):
        """Install needed dependencies from a list."""
        self.check_all()
        if deplist == ("all",):
            deplist = self.dependencies
        install_list = [
            dep
            for dep in deplist
            if not self.dependency_dict[dep]["installed"]
        ]
        if len(install_list):
            if not self.bin_path_exists:
                logger.error(
                    f"Installation directory {self.bin_path} does not"
                    " exist."
                )
                sys.exit(1)
            if not self.bin_path_writable:
                logger.error(
                    f"Installation directory {self.bin_path} is not"
                    " writable."
                )
                sys.exit(1)
            if not self.bin_path_in_path:
                logger.error(
                    f"Installation directory {self.bin_path} is not in"
                    " PATH."
                )
                sys.exit(1)
        for dep in install_list:
            self.install(dep)

    def _git(self, dep, dep_dict):
        """Git clone from list."""
        from sh import git  # isort:skip

        for repo in dep_dict["git_list"]:
            logger.debug(f"   cloning {dep} repo {repo}")
            git.clone(repo)

    def _download_untar(self, dep, dep_dict, progressbar=True):
        """Download and untar tarball."""
        from sh import tar  # isort:skip

        download_url = dep_dict["tarball"]
        dlname = download_url.split("/")[-1]
        download_path = Path(".") / dlname
        logger.debug(f"downloading {download_url} to {dlname}")
        tmp_path = download_path / (dlname + ".tmp")
        if progressbar:
            trackers = (ProgressTracker(DataTransferBar()),)
        else:
            trackers = None
        request_download(download_url, download_path, trackers=trackers)
        logger.debug(
            f"downloaded file {download_path}, size"
            f" {download_path.stat().st_size}"
        )
        tar_output = tar("xvf", download_path)
        logger.debug(tar_output)
        logger.debug("untar done")

    def _configure(self, dep, dep_dict):
        """Run make to build package."""
        logger.debug(f"   configuring {dep} in {Path.cwd()}")
        configure = sh.Command("./configure")
        try:
            configure_out = configure()
        except:
            logger.error("configure failed.")
            sys.exit(1)
        logger.debug(configure_out)

    def _copy_in_files(self, out_head, pkgname, depname, force=True):
        """Copy any files under installer_data to build directory."""
        resource_dir = "installer_data/" + depname
        if pkg_resources.resource_exists(pkgname, resource_dir):
            for root, dir, files in walk_package(pkgname, resource_dir):
                split_dir = os.path.split(root)
                if root == resource_dir:
                    out_subdir = ""
                else:
                    out_subdir = "/".join(list(split_dir)[2:])
                out_path = out_head / out_subdir
                if not out_path.exists() and len(files) > 0:
                    logger.info(f'Creating "{str(out_path)}" directory')
                    out_path.mkdir(0o755, parents=True)
                for filename in files:
                    data_string = pkgutil.get_data(
                        __name__, root + "/" + filename
                    ).decode("UTF-8")
                    file_path = out_path / filename
                    if file_path.exists() and not force:
                        logger.error(
                            f"File {str(file_path)} already exists."
                            + "  Use --force to overwrite."
                        )
                        sys.exit(1)
                    elif file_path.exists() and force:
                        operation = "Overwriting"
                    else:
                        operation = "Creating"
                    logger.debug(f"{operation} {file_path}")
                    with file_path.open(mode="wt") as fh:
                        fh.write(data_string)
                    if filename.endswith(".sh"):
                        file_path.chmod(0o755)

    def _make(self, dep, dep_dict):
        """Run make to build package."""
        from sh import make  # isort:skip

        logger.debug(f"   making {dep_dict['make']} in {Path.cwd()}")
        try:
            make_out = make(dep_dict["make"])
        except:
            logger.error("make failed.")
            logger.error(make_out)
            sys.exit(1)
        logger.debug(make_out)

    def _make_install(self, dep, dep_dict):
        """Run make install to install a package."""
        logger.info(f"   installing {dep} into {self.bin_path}")
        install_out = make.install(dep_dict["make_install"])
        logger.debug(install_out)

    def _copy_binaries(self, dep, dep_dict):
        """Copy the named binary to the bin directory."""
        logger.info(f"   copying {dep} into {self.bin_path}")
        for bin in dep_dict["copy_binaries"]:
            binpath = Path(bin)
            shutil.copy2(binpath, self.bin_path / binpath.name)

    def install(self, dep):
        """Install a particular dependency."""
        logger.info(f"installing {dep}")
        dep_dict = self.dependency_dict[dep]
        with tempfile.TemporaryDirectory() as tmp:
            tmppath = Path(tmp)
            logger.debug(f'build directory: "{tmppath}"')
            os.chdir(tmppath)
            #
            # Get the sources.  Alternatives are git or download
            #
            if "git_list" in dep_dict:
                self._git(dep, dep_dict)
            elif "tarball" in dep_dict:
                self._download_untar(dep, dep_dict)
            #
            # Change to the work directory.
            #
            logger.debug(f'building in directory {dep_dict["dir"]}')
            dirpath = Path(".") / dep_dict["dir"]
            if not dirpath.exists():
                raise ValueError(f'directory "{dirpath}" does not exist.')
            if not dirpath.is_dir():
                raise ValueError(f'directory "{dirpath}" is not a directory.')
            os.chdir(dirpath)
            workpath = Path.cwd()
            #
            # Copy in any additional files.
            #
            self._copy_in_files(workpath, self.pkg_name, dep, force=True)
            #
            # Build the executables.
            #
            if "configure" in dep_dict:
                self._configure(dep, dep_dict)
            if "configure_extra_dirs" in dep_dict:
                for newdir in dep_dict["configure_extra_dirs"]:
                    os.chdir(workpath / newdir)
                    self._configure(dep, dep_dict)
                    os.chdir(workpath)
            if "make" in dep_dict:
                self._make(dep, dep_dict)
            if "make_extra_dirs" in dep_dict:
                for newdir in dep_dict["make_extra_dirs"]:
                    os.chdir(workpath / newdir)
                    self._make(dep, dep_dict)
                    os.chdir(workpath)
            #
            # Install the executables.
            #
            if "make_install" in dep_dict:
                self._make_install(dep, dep_dict)
            elif "copy_binaries" in dep_dict:
                self._copy_binaries(dep, dep_dict)
