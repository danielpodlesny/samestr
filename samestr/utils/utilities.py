from datetime import datetime
from os.path import basename


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
