import re
import subprocess
from dataclasses import dataclass
from utils.exceptions import IpcressError

from adapters.primer3_to_ipcress import Primer3ToIpcressAdapter
from utils.write_output_files import write_ipcress_input
from utils.logger import Logger


@dataclass
class IpcressParams:
    pass


@dataclass
class IpcressResult:
    input_file: str
    stnd: bytes
    err: bytes


class Ipcress:
    def __init__(self, params) -> None:
        self.params = params
        self.logger = Logger(quiet=params["quiet"])

    def run(self) -> IpcressResult:
        self.logger.log("Running iPCRess...")

        self.logger.log('iPCRess params:')
        self.logger.log(self.params)
        params = self.params

        if params['primers']:                               # Comes from --primers argument, should be a path to a txt file
            self.logger.log('Loading custom iPCRess input file')
            input_path = params['primers']
            result = self.run_ipcress(input_path, params)
        else:
            self.logger.log('Building iPCRess input file.')

            adapter = Primer3ToIpcressAdapter()
            adapter.prepare_input(
                params['p3_csv'], params['min'], params['max'], params['dir']
            )
            input_path = write_ipcress_input(params['dir'], adapter.formatted_primers)

            result = self.run_ipcress(input_path, params)
            self.validate_primers(
                result.stnd.decode(), adapter.primer_data, params
            )

        self.logger.log('Finished!')

        return result

    def run_ipcress(self, input_path, params) -> IpcressResult:
        cmd = ' '.join([
            'ipcress', input_path, params['fasta'], '--mismatch', str(params['mismatch'])
        ])

        cmd = self.prettify_output(params['pretty'], cmd)

        self.logger.log(f'Running Exonerate iPCRess with the following command:\n{cmd}')

        result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stnd = result.stdout
        err = result.stderr
        if not stnd:
            raise IpcressError(err)

        return IpcressResult(input_path, stnd, err)

    @staticmethod
    def validate_primers(ipcress_output, primer_data, params, validate_coords=False) -> None:
        vp_logger = Logger(quiet=params["quiet"])
        if params["pretty"]:
            vp_logger.log('Output is pretty, skipping validation')
            return
        vp_logger.log('Validating primers...')

        for primer_pair in primer_data.keys():

            if validate_coords:
                fwd_coord = primer_data[primer_pair]['F']['start']
                rev_coord = primer_data[primer_pair]['R']['start']

                reg_exp = fr'{primer_pair} \d+ A {fwd_coord} 0 B {rev_coord} 0 forward'
            else:
                reg_exp = fr'{primer_pair} \d+ A \d+ 0 B \d+ 0 forward'

            match = re.search(reg_exp, ipcress_output)

            if not match:
                vp_logger.log(f'No valid primer pair found for {primer_pair}')

    @staticmethod
    def prettify_output(prettify, cmd):
        pretty_opt = 'false'
        if prettify:
            pretty_opt = 'true'

        return cmd + ' --pretty ' + pretty_opt
