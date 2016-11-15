import os
import yaml

def main(config_file='~/.SIDRE.yaml'):
    config_file = os.path.expanduser(config_file)
    if not os.path.exists(config_file):
        config = {'Default': True}
    else:
        with open(config_file, 'r') as FO:
            contents = FO.read()
            config = yaml.load(contents)
    print(config)


if __name__ == '__main__':
    main()