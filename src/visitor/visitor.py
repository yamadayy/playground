from .element import Directory, File
from abc import ABCMeta, abstractmethod
from functools import singledispatchmethod


class Visitor(metaclass=ABCMeta):
    @abstractmethod
    def visit(self, directory):
        pass


class ListVisitor(Visitor):
    def __init__(self):
        self.__currentdir = ''

    @singledispatchmethod
    def visit(self, directory):
        print("Error")

    @visit.register
    def _(self, directory:Directory):
        print("{0}/{1}".format(self.__currentdir, directory))
        savedir = self.__currentdir
        self.__currentdir = self.__currentdir + '/' + directory.getName()
        for f in directory:
            f.accept(self)
        self.__curentdir = savedir

    @visit.register
    def _(self, directory:File):
        print("{0}/{1}".format(self.__currentdir, directory))


"""
    def visit(self, directory):
        print("{0}/{1}".format(self.__currentdir, directory))
        if isinstance(directory, Directory):
            savedir = self.__currentdir
            self.__currentdir = self.__currentdir + '/' + directory.getName()
            for f in directory:
                f.accept(self)
            self.__curentdir = savedir
"""
