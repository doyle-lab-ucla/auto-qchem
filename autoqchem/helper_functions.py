import getpass
import glob
import logging
import os
from collections import Counter

import fabric
import paramiko

logger = logging.getLogger(__name__)


def ssh_connect(host, user) -> fabric.Connection:
    """Create ssh connection using fabric and paramiko, supports DUO authentication.

    :param host: remote host
    :param user: username to authenticate on remote host
    :return: fabric.Connection
    """

    client = paramiko.SSHClient()
    client.load_system_host_keys()

    def my_handler(title, instructions, prompt_list):
        print(f"{title}\n{instructions}")
        return [echo and input(prompt) or getpass.getpass(prompt) for (prompt, echo) in prompt_list]

    try:
        client.connect(host, username=user)
    except paramiko.ssh_exception.SSHException:
        pass

    client.get_transport().auth_interactive(username=user, handler=my_handler)

    c = fabric.Connection(host)
    c.client = client
    c.transport = client.get_transport()

    return c


def ssh_connect_password(host, user) -> fabric.Connection:
    """Create ssh connection using fabric and paramiko, password only (without DUO authentication).

    :param host: remote host
    :param user: username to authenticate on remote host
    :return: fabric.Connection
    """

    client = paramiko.SSHClient()
    client.load_system_host_keys()

    try:
        client.connect(host,username=user)
    except paramiko.ssh_exception.SSHException:
        pass

    client.get_transport().auth_password(username=user,password=getpass.getpass(f"{user}@{host}'s password:"))

    c = fabric.Connection(host)
    c.client = client
    c.transport = client.get_transport()

    return c


def ssh_connect_pem(host, user, pem_path) -> fabric.Connection:
    """Create ssh connection using fabric and paramiko, tries rsa key in .pem file, otherwise asks for password (without DUO authentication).

    :param host: remote host
    :param user: username to authenticate on remote host
    :param pem_path: private rsa key - full path to .pem file (optional)
    :return: fabric.Connection
    """

    client = paramiko.SSHClient()
    client.load_system_host_keys()

    try:
        if pem_path:
            client.connect(host,username=user,key_filename=pem_path)
        else:
            client.connect(host,username=user)
    except paramiko.ssh_exception.SSHException:
        pass

    if not pem_path:
        client.get_transport().auth_password(username=user,password=getpass.getpass(f"{user}@{host}'s password:"))
    
    c = fabric.Connection(host)
    c.client = client
    c.transport = client.get_transport()

    return c


def cleanup_directory_files(dir_path, types=()) -> None:
    """Remove files with specific extension(s) from a directory.

    :param dir_path: path of the directory to cleanup
    :param types: a tuple with file extenstions that will be removed
    """

    if not os.path.exists(dir_path):
        logger.debug(f"Directory {dir_path} does not exist. No files of type {types} to cleanup.")
        return

    count = 0
    for file_name in os.listdir(dir_path):
        if types:
            ext = os.path.splitext(file_name)[1][1:]
            if ext not in types:
                continue
        os.remove(os.path.join(dir_path, file_name))
        count += 1

    logger.debug(f"{count} files of type {types} removed from {dir_path}")


def cleanup_empty_dirs(dir_path) -> None:
    """Remove empty directories 1-level under the specified directory.

    :param dir_path: path of the directory to cleanup
    """

    for directory in glob.glob(f"{dir_path}/*/"):
        if not os.listdir(directory):
            os.rmdir(directory)


def convert_crlf_to_lf(file_path) -> None:
    """Convert windows line endings CRLF to unix line endings LF in a given file.

    :param file_path: file path to convert
    """

    windows_line_ending = b'\r\n'
    unix_line_ending = b'\n'

    # open file - need to use binary mode
    with open(file_path, 'rb') as open_file:
        content = open_file.read()

    # replace
    content = content.replace(windows_line_ending, unix_line_ending)

    # write
    with open(file_path, 'wb') as open_file:
        open_file.write(content)


def yes_or_no(question) -> bool:
    """Input question yes or no.

    :param question: question string
    :return: bool
    """

    while "the answer is invalid":
        reply = str(input(question + ' (y/n): ')).lower().strip()
        if reply[0] == 'y':
            return True
        if reply[0] == 'n':
            return False


def add_numbers_to_repeated_items(items_list) -> list:
    """Add numeric consecutive labels to repeated items in a list.

    :param items_list: list of strings
    :return: list
    """
    updated_items_list = []

    counts = Counter(items_list)
    current_count = {}
    for item in items_list:
        if counts[item] > 1:
            if item not in current_count:
                current_count[item] = 1
            updated_items_list.append(f"{item}{current_count[item]}")
            current_count[item] = current_count[item] + 1
        else:
            updated_items_list.append(item)
    return updated_items_list
