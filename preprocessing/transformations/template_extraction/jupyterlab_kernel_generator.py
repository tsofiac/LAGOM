#!/bin/env python3

import argparse
import grp
import os
import shutil
import stat
import subprocess
import sys
from datetime import datetime
from subprocess import PIPE, Popen

import jinja2
import yaml
from yaml import Loader

if sys.version_info[0] == 3 and sys.version_info[1] >= 8:
    from typing import TypedDict
else:
    from typing_extensions import TypedDict

ENCODING = "utf-8"

MINICONDA_ACTIVATE = "/projects/cc/se_users/carlsson_ksmq649/miniconda3/etc/profile.d/conda.sh"

CONDA_PREFIX = os.environ.get("CONDA_PREFIX") or None

JUPYTER_KERNELS = {
    "dev": "/opt/scp/services/jupyter/kernels-dev",
    "home": ".local",
    "prod": "/opt/scp/services/jupyter/kernels-prod",
}

REQUIRED_PACKAGES = ['ipykernel', 'jupyter']
USERS_ENVS = ["home"]
KERNEL_DIR = "share/jupyter/kernels"
BACKUP_DIR = "share/jupyter/kernels_backups"
SCRIPT_DIR = "/opt/scp/apps/system/software/jupyterlab_kernel_generator/1.0/" #os.path.dirname(os.path.abspath(__file__))
LUA_LMOD_PATH = os.path.join(
    os.getenv('MODULESHOME', '/usr/share/lmod/lmod'), 'libexec/lmod'
)
ADMIN_GROUPS = ["xem-scp-admins", "easybuild-jenkins"]
# code errors
ENOPKG = '65'  # Package not installed
EINVAL = '22'  # Invalid argument
EBADRQC = '56'  # Invalid request


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[32m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    OKMSG = ''.join([OKGREEN, "OK", ENDC])
    DELIM = ''.join([BOLD, '==========================', ENDC])


DOC_VENV = 'https://confluence.astrazeneca.com/x/HA_lDg#JupyterHUBrefresh-LocalKernels-Usingvirtualenv'
DOC_CONDA = 'https://confluence.astrazeneca.com/display/SCP/Using+conda'
DOC_SINGULARITY = 'https://confluence.astrazeneca.com/display/SCP/JupyterHUB+usage'


class KernelSpec(TypedDict):
    debug: bool
    display_name: str
    language: str
    manager: str
    packages: list


class JupyterKernel:
    def __init__(self, env, backup, display_name, manager, language, **kwargs):
        self.env = env
        self.backup = backup
        self.name = display_name
        self.manager = manager
        self.language = language
        self.kernel_path = self.get_kernel_path()
        self.kerneljson = os.path.join(self.kernel_path, "kernel.json")
        self.runsh = os.path.join(self.kernel_path, "run.sh")
        self.debug = kwargs.get('debug')
        self.kernel_map = self.create_kernel_map()

    def create_kernel_map(self):
        return {
            'language': self.language,
            'manager': self.manager,
            'display_name': self.name,
            'debug': self.debug,
            'kernel_path': self.runsh,
        }

    def create_kernel_dir(self) -> bool:
        """
        Creates kernel directory
        Return True if directory was created
        Return False if something went wrong
        """
        if os.path.exists(self.kernel_path) and self.backup:
            if not self.create_backup():
                return False
        try:
            os.umask(0o022)
            os.makedirs(self.kernel_path)
        except PermissionError:
            print(
                "Please, check that you have "
                f"write permissions at {self.kernel_path}"
            )
            return False
        except FileExistsError:
            print('Kernel directory exist')
            return False
        return True

    def create_backup(self) -> bool:
        """
        Moves files from Kernel Path to Backup Directory
        Return True if backup was created successfully
        Return False if something went wrong
        """
        current_time = str(datetime.now().strftime("%Y%m%d%H%M"))
        try:
            backup_path = f"{self.get_backup_path()}-{current_time}"
            shutil.move(self.kernel_path, backup_path)
            print(f"Backup was created: {backup_path}")
            return True
        except OSError:
            print(
                "You don't have permissions to move "
                f"{self.kernel_path} to {backup_path}"
            )
            return False

    def check_venv(self) -> bool:
        """
        Sanity check environment for kernel
        """
        return True

    def get_manager(self) -> str:
        """
        Return environment manager for the kernel
        """
        return self.manager

    def generate_json(self):
        """
        Generate kernel.json file
        """
        with open(os.path.join(SCRIPT_DIR, "templates/kernel.json")) as kernel_template:
            template = jinja2.Template(kernel_template.read())
            with open(self.kerneljson, "w") as json:
                json.write(template.render(**self.kernel_map))
        if not self.copy_icon():
            print("Failed to copy custom icon")
        print(f" {self.kerneljson} was generated")

    def generate_run(self):
        """
        Generate run.sh script

        """
        if self.debug == "":
            self.debug = False

        with open(os.path.join(SCRIPT_DIR, "templates/run.sh")) as run_template:
            template = jinja2.Template(run_template.read())
            with open(self.runsh, "w") as sh:
                sh.write(template.render(**self.kernel_map))

        self.fix_runsh_perms()
        print(f" {self.runsh} was generated")

    def get_kernel_path(self):
        """
        Return full path to kernel directory
        """
        if self.env == "home":
            return os.path.join(
                os.environ["HOME"], JUPYTER_KERNELS[self.env], KERNEL_DIR, self.name
            )

        return os.path.join(JUPYTER_KERNELS[self.env], KERNEL_DIR, self.name)

    def get_backup_path(self):
        """
        Return full path to kernel backup directory
        """
        if self.env == "home":
            return os.path.join(
                os.environ["HOME"], JUPYTER_KERNELS[self.env], BACKUP_DIR, self.name
            )
        return os.path.join(JUPYTER_KERNELS[self.env], BACKUP_DIR, self.name)

    def fix_runsh_perms(self) -> bool:
        """
        Set permissions to execute for user, group and other
        return True if add ugo+x to run.sh
        """
        try:
            st = os.stat(self.runsh)
            os.chmod(
                self.runsh, st.st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH
            )
        except OSError:
            return False

        return True

    def copy_icon(self) -> bool:
        """
        Copy custom icon to kernel path
        Return True if success
        """
        logo_64_path = os.path.join(self.kernel_path, "logo-64x64.png")
        try:
            shutil.copy(
                os.path.join(
                    SCRIPT_DIR, f"templates/icons/{self.language}/logo-64x64.png"
                ),
                logo_64_path,
            )
            # Fixed permissions
            os.chmod(logo_64_path, 0o444)
        except OSError:
            print(f"templates/icons/{self.language}/logo-64x64.png")
            return False
        print(f" {logo_64_path} copied")
        return True


class JupyterSingularity(JupyterKernel):
    def __init__(self, env, backup, display_name, manager, language, **kwargs):
        self.image_path = kwargs['image_path']
        super().__init__(env, backup, display_name, manager, language, **kwargs)

    def generate_run(self):

        self.kernel_map['image_path'] = self.image_path
        self.kernel_map['python_bin'] = 'python'

        super().generate_run()

    def check_venv(self) -> bool:
        if not os.path.exists(self.image_path):
            print(f'Singularity image {self.image_path} does not exists')
            return

        singularity_exec = f'singularity exec -e {self.image_path}'
        check_cmd = f'{singularity_exec} python -c "import ipykernel"'

        proc = subprocess.Popen(
            check_cmd,
            shell=True,
            executable="/bin/bash",
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        output, error = proc.communicate()
        if error:
            print('Can not find ipykernel module')
            print(f'Please, check documentation: {DOC_SINGULARITY}')
            return

        return True


class JupyterLmod(JupyterKernel):
    def __init__(self, env, backup, force, display_name, manager, language, **kwargs):
        super().__init__(env, backup, display_name, manager, language, **kwargs)
        self.force = force
        self.packages = kwargs.get('packages')
        self.debug = kwargs.get('debug')
        self.code_error = ''
        self.list_errors = {}

    def module(self, command, *arguments):
        numArgs = len(arguments)
        A = [LUA_LMOD_PATH, 'python', command]
        if numArgs == 1:
            A += arguments[0].split()
        else:
            A += list(arguments)

        proc = Popen(A, stdout=PIPE, stderr=PIPE)
        stdout, stderr = proc.communicate()
        stdout = stdout.decode("utf-8")
        stderr = stderr.decode("utf-8")
        err_out = sys.stderr
        self.code_error = proc.returncode
        exec(stdout)

    def check_venv(self):
        print(bcolors.DELIM)
        # load all modules from yml
        for module_name in self.packages:
            print('Loading lmod module -', module_name)
            output = JupyterLmod.module(self, 'load', module_name)
            message = bcolors.OKMSG
            if self.code_error != 0:
                message = ''.join(
                    [
                        bcolors.FAIL,
                        'WARNING: Cannot load module...\n',
                        module_name,
                        bcolors.ENDC,
                    ]
                )
                self.list_errors[module_name] = ENOPKG
            print(message)
        if len(self.list_errors) > 0 and not self.force:
            return False
        # checking if modules live
        print(bcolors.DELIM)
        for module_name in self.packages:
            print('Checking lmod module -', module_name)
            cmd = ''.join(
                [
                    'if ! module is-loaded ',
                    module_name,
                    '; then echo "WARNING: Module is not loaded..."; fi',
                ]
            )
            proc = subprocess.Popen(cmd, shell=True, stdout=PIPE)
            output = proc.communicate()[0].decode()
            message = bcolors.OKMSG
            if output:
                message = ''.join([bcolors.FAIL, output, module_name, bcolors.ENDC])
                self.list_errors[module_name] = EINVAL
            print(message)
        if len(self.list_errors) > 0 and not self.force:
            return False
        print(bcolors.DELIM)
        print(''.join([bcolors.OKBLUE, 'Checking ipykernel...', bcolors.ENDC]))
        # try to load ipykernel with new enviroment after loading modules
        cmd = 'SHORTPYV=`python -c \'import sys; print(".".join(map(str, sys.version_info[:2])))\'` && '
        cmd += 'PYTHONPATH=$EBROOTPYTHON/lib/python$SHORTPYV:$PYTHONPATH && python -c "import ipykernel"'

        try:
            ipykernelOut = subprocess.check_output(cmd, shell=True)
        except:
            self.list_errors['ipykernel'] = EBADRQC
            print(''.join([bcolors.FAIL, "ERROR: Critical, exit.", bcolors.ENDC]))
            return False
        print(bcolors.OKMSG)
        return True

    def generate_run(self):
        self.kernel_map['packages'] = self.packages
        self.kernel_map['python_bin'] = 'python'
        super().generate_run()


class JupyterConda(JupyterKernel):
    def __init__(self, env, backup, display_name, manager, language, **kwargs):
        super().__init__(env, backup, display_name, manager, language, **kwargs)
        self.venv = kwargs.get('conda_env')
        self.debug = kwargs.get('debug', False)
        self.kernel_map['python_bin'] = os.path.join(
            self.get_conda_prefix(), 'bin/python'
        )

        self.kernel_map['R_bin'] = os.path.join(self.get_conda_prefix(), 'lib/R/bin/R')

        self.kernel_map['kernel_cmd'] = 'IRkernel::main()'

    def show_conda_packages(self, packages):
        print(
            f"Find required packages.\nConda env '{self.venv}' include packages:\n"
            f"=================================="
        )
        for item in packages:
            print(f"{item}, version: {packages[item]} ")

    def get_conda_prefix(self):
        command = self.conda_activate_cmd()
        command.append(f"conda activate {self.venv}")
        command.append("echo $CONDA_PREFIX")
        cmd = " && ".join(command)

        proc = subprocess.Popen(
            cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE
        )
        output = proc.communicate()
        return output[0].decode(ENCODING).replace('\n', '')


    def show_conda_packages(self, packages):
        print(
            f"Find required packages.\nConda env '{self.venv}' include packages:\n"
            f"=================================="
        )
        for item in packages:
            print(f"{item}, version:{bcolors.OKBLUE} {packages[item]}{bcolors.ENDC}")
        print("==================================")

    def conda_activate_cmd(self):
        return [f"source {MINICONDA_ACTIVATE}"]

    def check_venv(self):
        print(f"Checking conda env '{self.venv}'...")

        command = self.conda_activate_cmd()
        command.append(f'conda activate {self.venv}')

        if self.language == 'python':
            command.append("conda list  'jupyter|ipykernel' ")

        if self.language == 'R':
            command.append("conda list  'IRkernel|jupyter'")
        cmd = " && ".join(command)
        proc = subprocess.Popen(
            cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE
        )
        output = proc.communicate()
        result = output[0].decode(ENCODING).split('\n')
        packages = {}
        for str in result:
            if not str.startswith("#") and str != '':
                name, version, *_ = str.split()
                packages[name] = version

        if self.language == 'python':
            if 'ipykernel' not in packages.keys() or 'jupyter' not in packages.keys():
                print(
                    "You need conda environment with ipykernel and jupyter packages.\n"
                    f"Check documentation: {DOC_CONDA}"
                )
                return False
            self.show_conda_packages(packages)
            return True

        if self.language == 'R':
            if 'jupyter' not in packages.keys() or 'r-irkernel' not in packages.keys():
                print(
                    "You need conda environment with r-irkernel and jupyter packages.\n"
                    f"Check documentation: {DOC_CONDA}"
                )
                return False
            self.show_conda_packages(packages)
            return True


class JupyterVenv(JupyterKernel):
    def __init__(self, env, backup, display_name, manager, language, **kwargs):
        super().__init__(env, backup, display_name, manager, language, **kwargs)
        self.venv_path = kwargs.get('venv_path')
        self.debug = kwargs.get('debug', False)

    def get_package_info(self, package) -> map:
        cmd = f"source {self.venv_path}/bin/activate && pip show {package}"
        proc = subprocess.Popen(
            cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE
        )
        output, _ = proc.communicate()
        result = output.decode(ENCODING)
        answer = {}
        for pkginfo in result.splitlines():
            try:
                key, value = pkginfo.split(": ")
            except ValueError:
                print(pkginfo.split(":"))
            answer[key] = value
        return answer

    def check_venv(self) -> bool:
        """
        Check virtual environment
        return False if incorrect venv path
        return False if didn't find package from REQUIRED_PACKAGES
        return True if find all packages from REQUIRED_PACKAGES
        """
        if os.path.exists(os.path.join(self.venv_path, 'bin/activate')):
            print(
                f"Find virual environment: {bcolors.OKBLUE} {self.venv_path} {bcolors.ENDC}"
            )
        else:
            print(f'{bcolors.FAIL} ERROR: Can not find environment {bcolors.ENDC}')
            return False

        for package in REQUIRED_PACKAGES:
            answer = self.get_package_info(package)
            if answer.get('Name') is not None:
                print(f'{answer["Name"]}: {answer["Version"]}')
            else:
                print(
                    f'Cannot find package: {package}.\nPlease, check documentation: {DOC_VENV}'
                )
                return False
        return True

    def generate_run(self):
        self.kernel_map['venv_path'] = self.venv_path
        self.kernel_map['python_bin'] = 'python'
        super().generate_run()


def is_admin() -> bool:
    """
    Check if user can be admin
    return TRUE is users group name is in ADMIN_GROUPS
    """
    gid = os.getgid()
    if grp.getgrgid(gid).gr_name in ADMIN_GROUPS:
        return True
    return False


def check_params(env: str) -> bool:
    """
    check if env parameter available
    return FALSE if env in unknown
    return FALSE if regular user try to use env not from USERS_ENVS
    """
    if env not in JUPYTER_KERNELS:
        kernels = list(JUPYTER_KERNELS.keys())
        print(f"Please, use one environment from list: {kernels}")
        return False
    if not is_admin() and env not in USERS_ENVS:
        print(f"You don't have access to environment: {env}")
        return False
    return True


if __name__ == "__main__":

    example_usage = """
    Examples:

    Generate kernel from definition file at home directory:
    jupyterlab_kernel_generator.py --definition config/Python-3.7.2-foss-2019a.yml

    Generate kernel from definition file at home directory without backup:
    jupyterlab_kernel_generator.py --definition config/Python-3.7.2-foss-2019a.yml --backup False
    """

    if is_admin():
        example_usage += """
    Generate kernel from definition file at dev directory:
    jupyterlab_kernel_generator.py --definition config/Python-3.7.2-foss-2019a.yml --environment dev
    """

    parser = argparse.ArgumentParser(
        description="Script for generating jupyterlab kernels",
        epilog=example_usage,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--environment",
        type=str,
        default="home",
        help=f"Provide environment destination (for {ADMIN_GROUPS} only)",
    )

    parser.add_argument(
        "--definition",
        metavar="def",
        type=argparse.FileType("rt"),
        help="Provide definition file",
        required=True,
    )

    parser.add_argument(
        "--backup", type=bool, default=True, help="Make backup before generate kernel"
    )

    parser.add_argument(
        "--force",
        default=False,
        action="store_true",
        help="Create kernel even if sanity-check doesn't pass",
    )

    # Parse arguments and set variables
    args = parser.parse_args()
    definition = args.definition
    env = args.environment
    backup = args.backup
    force = args.force

    try:
        kernelSpec = yaml.load(definition, Loader=Loader)
    except yaml.YAMLError as error:
        print(f"YAML error: {error}")
    except OSError as error:
        print(f"Could not open file: {error.strerror}")

    # Check environment
    if check_params(env):
        if kernelSpec["manager"] == "lmod":
            kernel = JupyterLmod(env, backup, force, **kernelSpec)
        elif kernelSpec["manager"] == "conda":
            kernel = JupyterConda(env, backup, **kernelSpec)
        elif kernelSpec["manager"] == "venv":
            kernel = JupyterVenv(env, backup, **kernelSpec)
        elif kernelSpec["manager"] == "singularity":
            kernel = JupyterSingularity(env, backup, **kernelSpec)
        else:
            print("Please, use lmod, conda or venv manager ")

        if kernel.check_venv():
            print("Checking environment: ", bcolors.OKMSG)
            if kernel.create_kernel_dir():
                kernel.generate_json()
                kernel.generate_run()
        else:
            print(f"Checking environment: {bcolors.FAIL} FAIL {bcolors.ENDC}")
