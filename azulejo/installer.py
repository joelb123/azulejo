# -*- coding: utf-8 -*-
"""Install binary dependencies."""
# standard library imports
import os
import pkg_resources
import pkgutil
import platform
import shutil
import sys
import tempfile
from pathlib import Path
from packaging import version

# third-party imports
import sh
from progressbar import DataTransferBar
from requests_download import ProgressTracker
from requests_download import download as request_download

# module imports
from .common import MinSpaceTracker
from .common import BUILD_DEV
from .common import free_mb
from .common import is_writable
from .common import logger

# Global constants
PLATFORM_LIST = ["linux", "bsd", "macos", "unknown"]
EXECUTABLE_EXTS = (".sh", ".pl", ".awk", ".py")
MEGABYTES = 1024.0 * 1024.0


def default_version_splitter(instring):
    """Split the version string out of version output."""
    return instring.split()[-1]


def walk_package(pkgname, root):
    """Walk through a package_resource."""
    dirs = []
    files = []
    for name in pkg_resources.resource_listdir(pkgname, str(root)):
        fullname = root / name
        if pkg_resources.resource_isdir(pkgname, str(fullname)):
            dirs.append(fullname)
        else:
            files.append(Path(name))
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
        accept_licenses=False,
        quiet=False,
    ):
        """Initialize dictionary of dependencies."""
        self.dependency_dict = dependency_dict
        self.force = force
        self.pkg_name = pkg_name
        self.dependencies = tuple(dependency_dict.keys())
        for dep in self.dependencies:
            if "required" not in self.dependency_dict[dep]:
                self.dependency_dict[dep]["required"] = True
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
        system = platform.system()
        if system == "Linux":
            self.platform = "linux"
        elif system.endswith("BSD"):
            self.platform = "bsd"
        elif system == "Darwin":
            self.platform = "macos"
        else:
            self.platform = "unknown"
            logger.warning(
                f"Unknown platform {system}, installs will likely fail"
            )
        self.not_my_platform = PLATFORM_LIST.copy()
        self.not_my_platform.remove(self.platform)
        self.accept_licenses = accept_licenses
        self.quiet = quiet
        if quiet:
            self.output_kwargs = {}
        else:
            self.output_kwargs = {"_out": sys.stdout, "_err": sys.stderr}
        self.status_msg = None

    def check_all(self, exe_paths=False):
        """Check all depenencies for existence and version."""
        self.status_msg = ""
        for dep in self.dependencies:
            target_version = version.parse(
                self.dependency_dict[dep]["version"]
            )
            version_command = self.dependency_dict[dep]["version_command"]
            self.dependency_dict[dep]["installed"] = not self.force
            for binary in self.dependency_dict[dep]["binaries"]:
                if sh.which(binary) == None:
                    self.dependency_dict[dep]["installed"] = False
                    exe = binary
                    ver_str = "not installed."
                else:
                    exe = sh.Command(binary)
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
                            f"at recommended version {installed_version}."
                        )
                    elif installed_version < target_version:
                        ver_str = (
                            f"installed {installed_version} <  target"
                            f" {target_version}."
                        )
                        self.dependency_dict[dep]["installed"] = False
                    elif installed_version > target_version:
                        ver_str = (
                            f"installed {installed_version} exceeds target"
                            f" {target_version}."
                        )
                if exe_paths:
                    self.status_msg += f"{exe} {ver_str}\n"
                else:
                    self.status_msg += f"{binary} {ver_str}\n"
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
        all_required = all(
            [
                self.dependency_dict[d]["installed"]
                for d in self.dependency_dict
                if self.dependency_dict[d]["required"]
            ]
        )
        all_optional = all(
            [
                self.dependency_dict[d]["installed"]
                for d in self.dependency_dict
                if not self.dependency_dict[d]["required"]
            ]
        )
        if all_required:
            if all_optional:
                optional_str = "and optional"
            else:
                optional_str = ""
            self.status_msg += (
                f"All required {optional_str} dependencies are installed.\n"
            )
        return all_required

    def status(self, exe_paths=False):
        """Returns and the installation status message."""
        if self.status_msg is None:
            self.check_all(exe_paths=exe_paths)
        return self.status_msg

    def install_list(self, deplist):
        """Install needed dependencies from a list."""
        self.check_all()
        if deplist == ("all",):
            deplist = [
                d
                for d in self.dependencies
                if self.dependency_dict[d]["required"]
            ]
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

    def _git(self, dep):
        """Git clone from list."""
        git = sh.Command("git")
        for repo in self.dependency_dict[dep]["git_list"]:
            logger.debug(f"   cloning {dep} repo {repo}")
            git.clone(repo, **self.output_kwargs)

    def _download_untar(self, dep):
        """Download and untar tarball."""
        tar = sh.Command("tar")
        download_url = self.dependency_dict[dep]["tarball"]
        dlname = download_url.split("/")[-1]
        download_path = Path(".") / dlname
        logger.debug(f"downloading {dep} at {download_url} to {dlname}")
        if self.quiet:
            trackers = ()
        else:
            trackers = (ProgressTracker(DataTransferBar()),)
        request_download(download_url, download_path, trackers=trackers)
        logger.debug(
            f"downloaded file {download_path}, size"
            f" {download_path.stat().st_size}"
        )
        try:
            tar("xvf", download_path, **self.output_kwargs)
        except:
            logger.error(f"untar of {download_path} failed")
            sys.exit(1)

    def _download_binaries(self, dep):
        """Download and untar tarball."""
        download_dict = self.dependency_dict[dep]["download_binaries"]
        if self.platform in download_dict:
            download_url = download_dict[self.platform]
        else:
            logger.warning(
                f"No binaries for download for {dep}, will fake it."
            )
            return
        dlname = download_url.split("/")[-1]
        download_path = Path(".") / dlname
        logger.debug(f"downloading {dep} at {download_url} to {dlname}")
        if self.quiet:
            trackers = ()
        else:
            trackers = (ProgressTracker(DataTransferBar()),)
        request_download(download_url, download_path, trackers=trackers)
        logger.debug(
            f"downloaded file {download_path}, size"
            f" {download_path.stat().st_size}"
        )

    def _configure(self, dep):
        """Run make to build package."""
        logger.debug(f"   configuring {dep} in {Path.cwd()}")
        configure = sh.Command("./configure")
        try:
            configure(**self.output_kwargs)
        except:
            logger.error(f"configure of {dep} failed.")
            sys.exit(1)

    def _copy_in_files(self, out_head, pkgname, dep, force=True):
        """Copy any files under installer_data to build directory."""
        resource_path = Path("installer_data") / dep
        platform_path = resource_path / self.platform
        not_my_platform_path_parts = [
            (resource_path / p).parts for p in self.not_my_platform
        ]
        if pkg_resources.resource_exists(pkgname, str(resource_path)):
            for root, unused_dirs, files in walk_package(
                pkgname, resource_path
            ):
                root_parts = root.parts
                if root_parts[:3] in not_my_platform_path_parts:
                    continue
                elif root_parts[:3] == platform_path.parts:
                    subdir_parts = root_parts[3:]
                else:
                    subdir_parts = root_parts[2:]
                out_path = out_head.joinpath(*subdir_parts)
                if not out_path.exists() and len(files) > 0:
                    logger.info(f'Creating "{str(out_path)}" directory')
                    out_path.mkdir(0o755, parents=True)
                for filename in files:
                    data_string = pkgutil.get_data(
                        __name__, str(root / filename)
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
                    if filename.suffix in EXECUTABLE_EXTS:
                        file_path.chmod(0o755)

    def _patch(self, dep):
        """Apply a patch."""
        patch = sh.Command("patch")
        logger.debug(
            f"   Applying patch {self.dependency_dict[dep]['patch']} in {Path.cwd()}"
        )
        try:
            patch(self.dependency_dict[dep]["patch"], **self.output_kwargs)
        except:
            logger.error(f"patch of {dep} failed.")
            sys.exit(1)

    def _make(self, dep):
        """Run make to build package."""
        make = sh.Command("make")
        logger.debug(
            f"   making {self.dependency_dict[dep]['make']} in {Path.cwd()}"
        )
        try:
            make(self.dependency_dict[dep]["make"], **self.output_kwargs)
        except:
            logger.error(f"make of {dep} failed.")
            sys.exit(1)

    def _make_install(self, dep):
        """Run make install to install a package."""
        make = sh.Command("make")
        logger.info(f"   installing {dep} into {self.bin_path}")
        try:
            make.install(
                self.dependency_dict[dep]["make_install"], **self.output_kwargs
            )
        except:
            logger.error(f"make install of {dep} failed.")

    def _copy_binaries(self, dep):
        """Copy the named binary to the bin directory."""
        logger.info(f"   copying {dep} into {self.bin_path}")
        for binary in self.dependency_dict[dep]["copy_binaries"]:
            binpath = Path(binary)
            shutil.copy2(binpath, self.bin_path / binpath.name)

    def _check_for_license_acceptance(self, dep):
        """Prompt for acceptance of license terms."""
        if "license" in self.dependency_dict[dep]:
            license_name = self.dependency_dict[dep]["license"]
        else:
            license_name = "restrictive"
        if "license_file" in self.dependency_dict[dep]:
            license_text = Path(
                self.dependency_dict[dep]["license_file"]
            ).read_text()
            logger.warning(license_text)
        while "invalid answer":
            reply = (
                str(
                    input(
                        f"Do you accept this {license_name} license? (y/n): "
                    )
                )
                .lower()
                .strip()
            )
            if len(reply) > 0:
                if reply[0] == "y":
                    return True
                if reply[0] == "n":
                    return False

    def install(self, dep):
        """Install a particular dependency."""
        logger.info(f"installing {dep}")
        dep_dict = self.dependency_dict[dep]
        build_dev_path = BUILD_DEV
        if not is_writable(build_dev_path):
            logger.error(f"Build directory {build_dev_path} is not writable")
            sys.exit(1)
        initial_free_mb = free_mb(build_dev_path)
        free_tracker = MinSpaceTracker(build_dev_path)
        logger.debug(
            f"building on {build_dev_path} with {free_tracker.report_min()} MB free"
        )
        with tempfile.TemporaryDirectory(
            prefix="installer-", dir=build_dev_path
        ) as tmp:
            tmppath = Path(tmp)
            logger.debug(f'build directory: "{tmppath}"')
            os.chdir(tmppath)
            #
            # Get the sources.  Alternatives are git or download
            #
            if "git_list" in dep_dict:
                self._git(dep)
            elif "tarball" in dep_dict:
                self._download_untar(dep)
            elif "download_binaries" in dep_dict:
                self._download_binaries(dep)
            free_tracker.check()
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
            # Check licensing
            #
            if "license" in dep_dict:
                if (
                    "license_restrictive" in dep_dict
                    and dep_dict["license_restrictive"]
                ):
                    logger.warning(
                        f"{dep} has a restrictive {dep_dict['license']} license"
                    )
                    if not self.accept_licenses:
                        accept = self._check_for_license_acceptance(dep)
                        if not accept:
                            logger.error(
                                f"License terms for {dep} not accepted, skipping installation"
                            )
                            return
                else:
                    logger.info(f"{dep} has a {dep_dict['license']} license")
            #
            # Build the executables.
            #
            if "patch" in dep_dict:
                self._patch(dep)
            if "configure" in dep_dict:
                self._configure(dep)
            if "configure_extra_dirs" in dep_dict:
                for newdir in dep_dict["configure_extra_dirs"]:
                    os.chdir(workpath / newdir)
                    self._configure(dep)
                    os.chdir(workpath)
            if "make" in dep_dict:
                self._make(dep)
            if "make_extra_dirs" in dep_dict:
                for newdir in dep_dict["make_extra_dirs"]:
                    os.chdir(workpath / newdir)
                    self._make(dep)
                    os.chdir(workpath)
            free_tracker.check()
            #
            # Install the executables.
            #
            if "make_install" in dep_dict:
                self._make_install(dep)
            elif "copy_binaries" in dep_dict:
                self._copy_binaries(dep)
            logger.info(
                f"Build of {dep} on {build_dev_path} used"
                + f" {free_tracker.report_used()}/{free_tracker.initial_space} MB"
            )
