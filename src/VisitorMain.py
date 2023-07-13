from visitor.visitor import ListVisitor
from visitor.element import File, Directory, FileTreatmentException


def start_main():
    try:
        print("makikng root entries")
        rootdir = Directory("root")
        bindir = Directory("bin")
        tmpdir = Directory("tmp")
        usrdir = Directory("usr")

        rootdir.add(bindir)
        rootdir.add(tmpdir)
        rootdir.add(usrdir)

        bindir.add(File("vi", 10000))
        bindir.add(File("latex", 20000))
        rootdir.accept(ListVisitor())

        print("")
        print("Making user entries...")
        yuki = Directory("yuki")
        hanako = Directory("hanako")
        tomura = Directory("tomura")

        usrdir.add(yuki)
        usrdir.add(hanako)
        usrdir.add(tomura)

        yuki.add(File("diary.html", 100))
        yuki.add(File("composite.py", 200))
        hanako.add(File("memo.tex", 300))
        tomura.add(File("game.doc", 400))
        tomura.add(File("junk.mail", 500))
        rootdir.accept(ListVisitor())

        print("")
        print("Occuring Exception..")
        tmpfile = File("tmp.txt", 300)
        bindir = Directory("bin")
        tmpfile.add(bindir)
    except FileTreatmentException as ex:
        print(ex.message)


if __name__ == '__main__':
    start_main()