from reader import get_args, read_config
from data_struct import Synthesis


def main(argv=''):
    args = get_args(argv)
    config = read_config(args.config)

    result = Synthesis(**config)
    print(result)


if __name__ == '__main__':
    main()