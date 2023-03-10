import csv
import json
from os import path, makedirs
from datetime import datetime
from typing import List

from utils.exceptions import FolderCreatorError, FileFormatError


class FolderCreator:
    dir = ''

    def set_dir(self, dir: str) -> None:
        self.dir = dir

    @classmethod
    def get_dir(self) -> str:
        return (self.dir)

    @classmethod
    def create(self, dir: str) -> None:
        try:
            if not path.isdir(dir):
                makedirs(dir)
                print(f'Folder {dir} is created')
            else:
                print(f'Warning: {dir} already exists and files may be overwritten')
        except OSError as error:
            raise FolderCreatorError(f'Unexpected OSError: {error}')

        self.set_dir(self, dir)

    @classmethod
    def create_timestamped(self, parent: str) -> None:
        PREFIX = 'td'

        dir_name = f'{PREFIX}_' + datetime.now().strftime('%Y%m%d%H%M%S%f')
        dir_full_path = path.join(parent, dir_name)

        self.create(dir_full_path)


def check_file_exists(file):
    if not path.exists(file):
        raise FileNotFoundError(f'Unable to find file: {file}')


def write_to_text_file(dir_path, data, file_name):
    file_path = dir_path + '/' + file_name + '.txt'

    file_h = {}
    if isinstance(data, list):
        file_h = open(file_path, "w")
        for row in data:
            file_h.write(row + "\n")
    else:
        file_h = open(file_path, "wb")
        file_h.write(data)
    file_h.close()

    print('Wrote to file: ' + file_path)

    return file_path


def read_csv_to_list_dict(csv_path, delimiter=',') -> List[dict]:
    check_file_exists(csv_path)

    data = []
    with open(csv_path, newline='') as csv_file:
        reader = csv.DictReader(csv_file, delimiter=delimiter)
        for row in reader:
            data.append(row)

    return data


def parse_json(file_path: str) -> dict:
    with open(file_path, "r") as file:
        try:
            result = json.load(file)
        except Exception as err:
            raise FileFormatError

    return result
