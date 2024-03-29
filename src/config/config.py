from abc import ABC, abstractmethod

from utils.file_system import parse_json


class Config(ABC):
    def __init__(self):
        pass

    @abstractmethod
    def read_config(self):
        pass
        

class DesignerConfig(Config):
    def __init__(self, config_file: str = None):
        self._default_config_file = 'config/designer.config.json'

        config = self.read_config(self._default_config_file, config_file)
        self.params = {'stringency_vector': config['stringency_vector']}

    @staticmethod
    def read_config(
            default_config_file: str,
            config_file: str = None,
    ) -> dict:
        default_config = parse_json(default_config_file)
        keys = default_config.keys()

        if config_file is None or config_file == default_config_file:
            return default_config
        else:
            config = parse_json(config_file)

            for field in keys:
                config.setdefault(field, default_config[field])

        return config


class Primer3ParamsConfig(Config):
    def __init__(self, config_file: str = None):
        self._default_config_file =  'src/primer/primer3.config.json'

        self.params = self.read_config(self._default_config_file, config_file)

    @staticmethod
    def read_config(
            default_config_file: str,
            config_file: str = None,
    ) -> dict:

        file = default_config_file
        if config_file is not None:
            file = config_file

        config_data = parse_json(file)

        return config_data
