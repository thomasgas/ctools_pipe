from environs import Env
import os


def create_path(path_to_expand):
    env = Env()
    out_path = path_to_expand

    # try just to interpret the env variable (removing the dollar)
    # if it's $VAR/folder then split
    if out_path.startswith("$"):
        out_path = out_path[1:]
        try:
            out_path = env(out_path)
        except:
            env_var = out_path.split('/', 1)[0]
            folders = out_path.split('/', 1)[1]
            out_path = env(env_var) + '/' + folders

    if out_path.endswith("/"):  # just used to remove backslash if applied
        out_path = out_path[:-1]

    if not os.path.exists(out_path):
        os.makedirs(out_path)

    return out_path
