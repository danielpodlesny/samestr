import os
import re

from datetime import datetime
from os.path import basename, isdir, isfile
import numpy as np
import pandas as pd
import gzip
import json


def log_time(method):
    """
        Decorator for timing of function execution.
    """
    def timer(*args, **kw):
        t_start = datetime.now()
        result = method(*args, **kw)
        t_end = datetime.now()
        s = int((t_end - t_start).seconds)
        hours, remainder = divmod(s, 3600)
        minutes, seconds = divmod(remainder, 60)
        t_delta = '{:02}h:{:02}m:{:02}s'.format(int(hours), int(minutes),
                                                int(seconds))
        if 'log_time' in kw:
            name = kw.get('log_name', method.__name__.upper())
            kw['log_time'][name] = t_delta
        else:
            if minutes < 1:
                print(('Finished "{}" in {:2}s'.format(method.__name__,
                                                       int(seconds))))
            elif hours < 1:
                print(('Finished "{}" in {:2}m {:2}s'.format(
                    method.__name__, int(minutes), int(seconds))))
            else:
                print(('Finished "{}" in {:2}h {:2}m {:2}s'.format(
                    method.__name__, int(hours), int(minutes), int(seconds))))

        return result

    return timer


def is_exe(name):
    """
        Check whether `name` is on PATH.
    """

    from distutils.spawn import find_executable
    return find_executable(name) is not None


def all_exe(exe_list):
    """
        Check whether each `exe` in `exe_list` is on PATH.
    """

    for exe in exe_list:
        if not is_exe(exe):
            raise Exception('Executable for %s not found' % exe)
        print(('Found in path: %s') % exe)

    return True


def list_str_all_endswith(string_list, ending):
    """
        Returns True if all `list` items share the same `ending`
    """
    if not all([item.endswith(ending) for item in string_list]):
        return False
    return True


def str_cut_endswith(string, ending):
    """
        Cuts `string` at the `ending` position
    """
    if string.endswith(ending):
        return string[:-len(ending)]
    return string


def list_str_cut_endswith(string_list, endings):
    """
        Returns list of strings with `endings` removed
    """
    if type(endings) is not list:
        endings = list(endings)

    cut_strings = []
    for s in string_list:
        cut = False
        for e in endings:
            if s.endswith(e):
                cut_strings.append(str_cut_endswith(s, e))
                cut = True
        if not cut:
            cut_strings.append(s)

    return cut_strings


def list_group_by_basename(names, cut_name_endings):
    cut_names = list_str_cut_endswith(names, cut_name_endings)
    base_names = {}
    for cut_name, input_file in zip(cut_names, names):
        cut_name = basename(cut_name)
        if cut_name in base_names:
            base_names[cut_name].append(input_file)
        else:
            base_names[cut_name] = [input_file]
    return base_names

def load_numpy_file(input_file):
    """
        Load numpy file
    """
    if input_file.endswith('.npz'):
        return np.load(input_file, allow_pickle=True)['arr_0']
    elif input_file.endswith('.npy.gz'):
        with gzip.open(input_file, 'rb') as f:
            return np.load(f, allow_pickle=True)
    else:
        return np.load(input_file, allow_pickle=True)

def collect_input_files(input_dir, accepted_extensions, recursive=False):
    """
    Collect filepaths with valid extensions from an input directory.
    If recursive is True, collect files from subdirectories as well.
    """

    if not os.path.exists(input_dir):
        raise ValueError(f"Could not find input directory {input_dir} ({os.path.abspath(input_dir)})")

    samples = {}
    extensions = set()
    accepted_extensions = {ext.strip(".") for ext in accepted_extensions}

    for sample_dir, _, files in os.walk(input_dir):
         
        for fn in files:
            extension = None

            # accepted extensions can be single ('.fa') or composite ('.fa.gz')
            file_extensions = re.split(r"\.", fn)[-2:]
            if len(file_extensions) >= 1:
                # -> check last or last two suffixes                
                if file_extensions[-1] in accepted_extensions:
                    # check rightmost suffix
                    extension = file_extensions[-1]
                elif len(file_extensions) == 2:
                    # else check concatenation of two rightmost suffixes
                    ext = ".".join(file_extensions)
                    if ext in accepted_extensions:
                        extension = ext

            if extension is not None:
                # if this is a file with a valid extension
                extension = f".{extension}"
                extensions.add(extension)
                if len(extensions) > 1:
                    # we only allow one extension per input set
                    # error out if more than one found
                    raise ValueError(f"Found files with different extensions: {extensions}")
                # strip extension off filename to obtain sample id
                sample = fn[:-len(extension)]
                # collect file
                samples.setdefault(sample, []).append(os.path.join(sample_dir, fn))

        if not recursive:
            break

    return samples, (list(extensions) + [None])[0]

def read_json(fpath):
    with open(fpath, 'r') as file:
        return json.load(file)
    
def write_json(fpath, data):
    with open(fpath, 'w') as json_file:
        json.dump(data, json_file, indent=4)

def read_tsv(file_path):
    return pd.read_csv(file_path, sep='\t')

def write_tsv(file_path, df):
    df.to_csv(file_path, sep='\t', index=False)

def load_samestr_db_manifest(db_path):
    """
    read the samestr database files:
    - db_manifest.json
    - db_clades.json
    - db_taxonomy.tsv
    """

    db_manifest_fn = db_path + '/' + 'db_manifest.json'
    db_clades_fn = db_path + '/' + 'db_clades.json'
    db_taxonomy_fn = db_path + '/' + 'db_taxonomy.tsv'

    if isfile(db_manifest_fn):
        db_existed = True
        db_manifest = read_json(db_manifest_fn)
        db_clades = read_json(db_clades_fn)
        db_taxonomy = {'fpath': db_taxonomy_fn,
                       'records': read_tsv(db_taxonomy_fn)}
    else:
        db_existed = False
        db_manifest = {'fpath': db_manifest_fn,
                       'database':{},
                       'records':{'total_n_files':0, 'total_n_clades':0, 'total_n_markers': 0, 'total_n_positions': 0}}
        db_clades = {'fpath': db_clades_fn,
                     'records':{}}
        db_taxonomy = {'fpath': db_taxonomy_fn,
                       'records': pd.DataFrame(columns=['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'clade']) }
        
    # manifest existed, manifest
    return db_existed, db_manifest, db_clades, db_taxonomy
