import re
from argparse import ArgumentParser
from configparser import ConfigParser
from os import path


def get_args(argv):
    parser = ArgumentParser()
    parser.add_argument('config', nargs='?', default='config.ini')
    if argv:
        if isinstance(argv, str):
            argv = argv.split()
        elif not hasattr(argv, '__iter__'):
            raise TypeError(f'argv must be `str` or iterable, not {type(argv)}')
        args = parser.parse_args(argv)
    else:
        # get arguments from terminal
        args = parser.parse_args()
    assert path.isfile(args.config), f'config: `{args.config}` must be an existing file'
    return args


def parse_line(line):
    words = re.split(' *, *| *: +| +', line)
    key, *vals = words
    for i, v in enumerate(vals):
        if re.match('[-+]?\d+$', v):
            v = int(v)
        elif re.match('(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?$', v):
            v = float(v)
        elif v.lower().strip() in ['yes', 'no', 'on', 'off', 'true', 'false', '1', '0']:
            v = v.lower().strip() in ('yes', 'on', 'true', '1')
        vals[i] = v
    if vals:
        if len(vals) == 1:
            value = vals[0]
        else:
            value = vals
    else:
        value = ''
    return key, value


def parse_config(config):
    # TODO: handle errors to reflect the name of the field that failed
    # TODO: create functions to handle reading sections
    config_dict = dict()
    for name, value in config['CATALYST'].items():
        if ' ' in name:
            name = name.replace(' ', '_')
        if name == 'percentage':
            if value[0] == 'a':
                res = 'atomic'
            elif value[0] == 'w':
                res = 'weight'
            else:
                raise ValueError(f"percentage can only be 'a' or 'w', not {value}")
        elif name == 'support':
            support, percentage = value.split()
            res = (support, float(percentage))
        elif '\n' in value:
            res = dict()
            for line in value.strip().split('\n'):
                k, v = parse_line(line)
                res[k] = v
        else:
            res = float(value)
        config_dict[name] = res

    config_dict['experimental'] = dict()
    experimental_config = config._sections.get('EXPERIMENTAL', dict())
    for name, value in experimental_config.items():
        if name == 'precursors':
            res = dict()
            for line in value.strip().split('\n'):
                k, v = parse_line(line)
                res[k] = v
        elif name == 'catalyst_mass':
            res = float(value)
        elif name == 'support':
            support, percentage = value.split()
            res = (support, float(percentage))
        config_dict['experimental'][name] = res

    ink = dict()
    ink_config = config._sections.get('INK', dict())
    for name, value in ink_config.items():
        if name == 'PGM':
            v = value.split()
        else:
            v = float(value)
        ink[name] = v
    config_dict['ink'] = ink

    return config_dict


def read_config(config_file):
    config = ConfigParser(interpolation=None)
    config.read(config_file)

    config = parse_config(config)
    return config
