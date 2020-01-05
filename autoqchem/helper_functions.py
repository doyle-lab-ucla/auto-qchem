import getpass
import logging
import os

import fabric
import paramiko

logger = logging.getLogger(__name__)


def cleanup_directory_files(dir_path, types=()) -> None:
    """remove files of specific extension(s) from a directory"""

    if not os.path.exists(dir_path):
        logger.info(f"Directory {dir_path} does not exist. No files of type {types} to cleanup.")
        return

    count = 0
    for file_name in os.listdir(dir_path):
        if types:
            ext = os.path.splitext(file_name)[1][1:]
            if ext not in types:
                continue
        os.remove(os.path.join(dir_path, file_name))
        count += 1

    logger.info(f"{count} files of type {types} removed from {dir_path}")


def convert_crlf_to_lf(file_path) -> None:
    """utility to convert a windows line endings CRLF (\r\n)
    to unix line endings LF (\n)"""

    # replacement strings
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
    """input question yes or no"""

    while "the answer is invalid":
        reply = str(input(question + ' (y/n): ')).lower().strip()
        if reply[0] == 'y':
            return True
        if reply[0] == 'n':
            return False


def ssh_connect(host, user) -> fabric.Connection:
    """ssh connection using fabric and paramiko
    special handling is required to manage DUO
    authentication"""

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

    c = fabric.Connection("")
    c.client = client
    c.transport = client.get_transport()

    return c