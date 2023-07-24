from astropy.constants import Constant
import pkg_resources
import yaml
import codecs
import numpy as np
import importlib


class _TemporaryConstant(Constant):
    # TODO 規格の切り返したい場合は、このクラスから切り離した方が良いです
    # 規格の切り替えはこの辺で行ってます
    # see https://onl.sc/BTGXenG (astropyのconfig.py)
    default_reference = "Temporary Constant"
    _registry = {}
    _has_incompatible_units = set()


_formulas = set([])


def constant_formula(func):
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)
    _formulas.add(wrapper)
    return wrapper


class Parameters:
    __instance = None
    __ignore_list = (
        '__is_dirty',
        '__check_dirty',
        '__constants',
        '__load_file',
        'ready',
        'apply',
        '_translate_value',
        '_extract_value'
    )

    def __new__(cls, *args, **kwargs):
        if Parameters.__instance is None:
            return super(Parameters, cls).__new__(cls)
        return Parameters.__instance

    def __init__(self):
        if Parameters.__instance is not None:
            return
        self.__is_dirty = False  # TODO ファイル読み込みしないならここでTrueにして良い
        self.__constants = {}
        # TODO 定数の規格を切り替えるならimportlibを一工夫する
        const = importlib.import_module("jasmine_toolkit.utils.constants.jasmine_constant")
        # print(dir(const))
        for name in dir(const):
            if name.startswith("__"):
                continue
            self.__constants[name] = getattr(const, name)
        # filename = pkg_resources.resource_filename(
        #     'jasmine_toolkit', 'utils/constants/constants.yaml')
        # self.__load_file(filename, True)
        # print(self.__constants)
        self.__is_dirty = True
        Parameters.__instance = self

    def __setattr__(self, name, value):
        if not name.endswith(Parameters.__ignore_list):
            self.__check_dirty(False)
            if name in self.__constants:
                self.__constants[name] = value
                return
        object.__setattr__(self, name, value)

    def __getattr__(self, name):
        # print(name)
        if not name.endswith(Parameters.__ignore_list):
            self.__check_dirty(True)
            if name in self.__constants:
                return self.__get_value(name)
        return object.__getattribute__(self, name)

    def __get_value(self, name):
        ret = self.__constants[name]
        try:
            if ret not in _formulas:
                return ret
        except TypeError:
            return ret
        return ret()

    def is_dirty(self):
        return self.__is_dirty

    def ready(self):
        self.__is_dirty = False

    def apply(self, filename):
        self.__load_file(filename, False)

    def __load_file(self, filename, init):
        with codecs.open(filename, encoding='utf-8') as file:
            obj = yaml.safe_load(file)
            # print(obj)
            for name, val in obj.items():
                if init or name in self.__constants:
                    v = self._translate_value(name, val)
                    self.__constants[name] = v
                else:
                    cls_name = __class__.__name__
                    raise AttributeError(
                        f"'{cls_name}' object has no attribute '{name}'")

    def __check_dirty(self, check):
        if self.__is_dirty == check:
            if check:
                raise RuntimeError('getter don\'t call!!')
            else:
                raise RuntimeError('setter don\'t call!!')

    def _translate_value(self, key, dic):
        description = _extract_description(dic)
        description = description if not description == '' else key
        unit = _extract_unit(dic)
        value = self._extract_value(dic)
        if type(value) == float or type(value) == int:
            # TODO uncertainty support.
            return _TemporaryConstant(key, description, value, unit, 0.0)
        return value

    def _extract_value(self, dic):
        val = _extract_value(dic, 'value')
        if not type(val) is str:
            return val
        try:
            return float(val) if _maybe_real(val) else int(val)
        except ValueError:
            pass
        try:
            # print(val)
            v = eval(val)
            # print(v)
            # print(type(v))
            return v
        except BaseException as e:
            # print(e)
            return val


def _maybe_real(val):
    return '.' in val or 'e' in val.lower()


def _extract_description(dic):
    return _extract_value(dic, 'description')


def _extract_unit(dic):
    return _extract_value(dic, 'unit')


def _extract_value(dic, key):
    ret = dic[key] if key in dic else ''
    return ret if ret is not None else ''


Parameters()   # Python 3.7以降なら__new__自体がスレッドセーフなのでいらないらしいです
