import sys
import os
import ConfigParser

def initialize_settings(config_file, project=None):
    """Stores config file info into dictionary"""
    config = ConfigParser.ConfigParser()
    config.optionxform=str # preserve uppercase
    config.read(config_file)

    settings = {}             
    for section in config.sections():
        settings[section] = {}
        for name in config.options(section):
            final_section = section
            value = config.get(section, name)
            
            if project is not None:
                if project == section:
                    if name.endswith('-cmd'):
                        final_section = 'commands'
                        name = name[:-4]
                
                    elif name.endswith('-mem'):
                        final_section = 'memory'
                        name = name[:-4]
                    
                    elif name.endswith('-tmpmem'):
                        final_section = 'tmpmem'
                        name = name[:-7]

                    elif name.endswith('-java'):
                        final_section = 'java'
                        name = name[:-5]
            
            settings[final_section][name] = value

    return settings

def get_value(values, section, name):
    """Given dictionary(values), section, and name, returns value"""
    if values.has_key(section) and values[section].has_key(name):
        return values[section][name]
    else:
        return None
    
def create_args(values, section, name, args_dict=None):
    """Replaces variables in config file argments with values in dictionary"""
    args = get_value(values, section, name)    
    if args:
        for field,value in args_dict.iteritems():
            args = args.replace('${' + str(field) + '}', str(value))
    else:
        args = ''
            
    return args

