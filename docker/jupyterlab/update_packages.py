import subprocess
import requests
import toml
import re
import importlib.metadata


def extract_git_url(line):
    """Extracts a git URL and optional branch name from a line."""
    print(line)
    match = re.search(r"https://github\.com/([^/@]+/[^/@]+)(?:@([^#]+))?", line)
    print(match)
    print(match.groups())
    if match:
        repo, branch = match.groups()
        branch = branch or "main"  # Default to 'main' if branch is not specified
        print(repo, branch)
        return f"https://github.com/{repo}", branch

def add_prefix_to_numbers(input_string):
    # Regular expression to find numbers
    pattern = r'\d+'
    
    # Function to add "%23" before each number
    def add_prefix(match):
        return f"%23{match.group(0)}"
    
    # Replace each number in the string with "%23" followed by the number
    modified_string = re.sub(pattern, add_prefix, input_string)
    
    return modified_string

def get_pyproject_toml_url(git_url, branch="main"):
    """Constructs the URL to the raw pyproject.toml file."""
    repo_path = git_url.replace("https://github.com/", "")
    repo_path = repo_path.replace(".git", "")
    # careful might brake if the branch name is not issue
    # check if a number is in the branch name, and add before the number present %23
    branch = add_prefix_to_numbers(branch)
        
    return f"https://raw.githubusercontent.com/{repo_path}/{branch}/pyproject.toml"


def compare_installed_packages(package_name, version):
    """Compares listed packages and their versions with the ones installed in the current environment."""
    try:
        installed_version = importlib.metadata.version(package_name)
        if installed_version == version:
            print(
                f"{package_name} local=={version} remote=={installed_version} - No update needed"
            )
            return False

        else:
            print(
                f"{package_name} local=={version} remote=={installed_version} - Update needed"
            )
            return True
    except importlib.metadata.PackageNotFoundError:
        print(f"{package_name} is not installed")
        return True


def get_version_from_git_url(line):
    """Extracts and returns the version number from a pyproject.toml given a line with a GitHub repository URL."""
    git_url, branch = extract_git_url(line)
    if git_url:
        pyproject_toml_url = get_pyproject_toml_url(git_url, branch)
        print("url is:", pyproject_toml_url)
        response = requests.get(pyproject_toml_url)
        if response.status_code == 200:
            try:
                pyproject_content = toml.loads(response.text)
                package_name = pyproject_content["tool"]["poetry"]["name"]
                version = pyproject_content["tool"]["poetry"]["version"]
                print(f"Found version {version} for {package_name} in pyproject.toml")
                return package_name, version, git_url, branch
            except KeyError:
                print(f"Version information not found in pyproject.toml for {git_url}")
        else:
            print(
                f"Failed to fetch pyproject.toml from {git_url}: HTTP {response.status_code}"
            )
    else:
        print(f"Invalid git URL: {line}")


def install_package(package_name, version=None, git_url=None, branch="main"):
    # Determine if it's a pip package or a git repository
    is_pip = git_url is None

    install_command = ["pip", "install", "-q", "--force-reinstall"]

    if is_pip:
        # For PyPI packages
        print(f"Installing {package_name}=={version} from PyPI...")
        package_spec = f"{package_name}=={version}" if version else package_name
        install_command.append(package_spec)
    else:
        # For Git repositories
        print(
            f"Installing {package_name}=={version} from {git_url} on branch {branch} ..."
        )
        branch_spec = f"@{branch}" if branch else ""
        if ".git" not in git_url:
            git_url = f"{git_url}.git"
        git_install_url = f"git+{git_url}{branch_spec}"
        print('git_install_url:', git_install_url)
        install_command.append(git_install_url)

    # Execute the install command
    subprocess.run(install_command, check=True)


def process_requirements_file(requirements_path):
    print(f"checking for updates in {requirements_path}")
    """Processes a requirements.txt file to print version numbers for git repositories."""
    with open(requirements_path, "r") as file:
        for line in file:
            print('line is:', line)
            if not line.strip():
                continue
            line = line.strip()
            if "==" in line:
                package_name, version = line.split("==")
                git_url = None
                branch = None
            else:
                package_name, version, git_url, branch = get_version_from_git_url(line)
                print('package_name:', package_name, 'version:', version, 'git_url:', git_url, 'branch:', branch)

            reinstall = compare_installed_packages(package_name, version)

            if reinstall:
                install_package(package_name, version, git_url, branch=branch)


if __name__ == "__main__":
    process_requirements_file("/tmp/requirements.txt")
