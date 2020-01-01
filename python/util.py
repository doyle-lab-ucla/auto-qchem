import os
import logging

logger = logging.getLogger(__name__)


def cleanup_directory_files(dir_path, types=[]):
    """remove files of sepcific extension(s) from a directory"""

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


def convert_crlf_to_lf(file_path):
    """utilitiy to convert a windows line endings CRLF (\r\n)
    to unix line endings LF (\n)"""

    # replacement strings
    WINDOWS_LINE_ENDING = b'\r\n'
    UNIX_LINE_ENDING = b'\n'

    # open file - need to use binary mode
    with open(file_path, 'rb') as open_file:
        content = open_file.read()

    # replace
    content = content.replace(WINDOWS_LINE_ENDING, UNIX_LINE_ENDING)

    # write
    with open(file_path, 'wb') as open_file:
        open_file.write(content)

def yes_or_no(question):
    """input question yes or no"""

    while "the answer is invalid":
        reply = str(input(question+' (y/n): ')).lower().strip()
        if reply[0] == 'y':
            return True
        if reply[0] == 'n':
            return False