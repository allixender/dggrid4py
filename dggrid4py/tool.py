import requests
import stat
import os

def download_executable(url, folder="./"):
    # Create the directory if it doesn't exist
    os.makedirs(folder, exist_ok=True)
    
    # Get the filename from the URL
    filename = url.split('/')[-1]
    local_path = os.path.join(folder, filename)
    
    # Download the file
    response = requests.get(url, stream=True)
    response.raise_for_status()  # Raise an exception for HTTP errors
    
    # Save the file
    with open(local_path, 'wb') as f:
        for chunk in response.iter_content(chunk_size=8192):
            f.write(chunk)
    
    # Make the file executable (add execute permission)
    current_permissions = os.stat(local_path).st_mode
    os.chmod(local_path, current_permissions | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    
    # Return the absolute path
    return os.path.abspath(local_path)

# TODO: could be cache checked
def get_portable_executable(folder="./"):
    import platform

    uname = platform.uname()
    cpu = uname.machine.lower()
    opsystem = uname.system.lower()
    
    files = {
        'linux-aarch64': 'https://github.com/allixender/dggrid4py/releases/download/v0.3.2/dggrid-linux-aarch-gnu',
        'linux-x86_64': 'https://github.com/allixender/dggrid4py/releases/download/v0.3.2/dggrid-linux-x64',
        'darwin-arm64': 'https://github.com/allixender/dggrid4py/releases/download/v0.3.2/dggrid-macos-aarch',
        'windows-amd64': 'https://github.com/allixender/dggrid4py/releases/download/v0.3.2/dggrid-windows-x64.exe'
    }

    url = files[f"{opsystem}-{cpu}"]
    if url is None:
        raise ValueError(f"No portable executable available for {opsystem} {cpu}. Please use dggrid4py with a local DGGRID installation.")
    file_path = download_executable(url, folder)
    return file_path
