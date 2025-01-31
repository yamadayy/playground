import numpy as np


class Star:
    def __init__(self, x, id):
        self.position = x
        self.id = id
        self.next = None

    def set_id(self, id):
        self.id = id


class Stars:
    def __init__(self):
        self.head = Star(None, None)

    def insert(self, x, id):
        front = self.head
        rear = front.next
        while rear and id > rear.id:
            front = rear
            rear = rear.next
        star = Star(x, id)
        star.next = rear
        front.next = star

    def delete(self, id):
        front = self.head
        rear = front.next
        while rear:
            if rear.id == id:
                break
            front = rear
            rear = rear.next
        if not rear:
            print("[*] Data not found")
            return
        front.next = rear.next
        rear = None

    def reset_id(self):
        tmp = self.head.next
        idx = 0
        while tmp:
            tmp.set_id(idx)
            tmp = tmp.next
            idx = idx + 1

    def show(self):
        tmp = self.head.next
        while tmp:
            print("{},{}".format(tmp.id, tmp.position))
            tmp = tmp.next


class Plate:
    def __init__(self, xmin, xmax, id):
        self.xmin = xmin
        self.xmax = xmax
        self.id = id
        self.next = None

    def set_id(self, id):
        self.id = id


class Plates:
    def __init__(self):
        self.head = Plate(None, None, None)

    def insert(self, xmin, xmax, id):
        front = self.head
        rear = front.next
        while rear and id > rear.id:
            front = rear
            rear = rear.next
        star = Plate(xmin, xmax, id)
        star.next = rear
        front.next = star

    def delete(self, id):
        front = self.head
        rear = front.next
        while rear:
            if rear.id == id:
                break
            front = rear
            rear = rear.next
        if not rear:
            print("[*] Data not found")
            return
        front.next = rear.next
        rear = None

    def reset_id(self):
        tmp = self.head.next
        idx = 0
        while tmp:
            tmp.set_id(idx)
            tmp = tmp.next
            idx = idx + 1

    def show(self):
        tmp = self.head.next
        while tmp:
            print("{},{}~{}".format(tmp.id, tmp.xmin, tmp.xmax))
            tmp = tmp.next


def make_stars(num_stars):
    tmp = np.zeros(num_stars)
    for l in range(num_stars):
        tmp[l] = 0.5 + (l + 0.5) / 2 / num_stars
    return tmp


def make_plate():
    tmp = np.zeros(2)
    tmp[0] = 0.5
    tmp[1] = 1.0
    return tmp


def make_design_matrix(num_stars, plate_dependent_max_power, plate_independent_max_power, stars, plates):
    tmp = np.zeros((2 * num_stars, num_stars + plate_dependent_max_power + plate_independent_max_power + 2))
    for l in range(num_stars):
        tmp[l][l] = 1
        tmp[l + num_stars][l] = 1
        tmp[l][num_stars] = -1
        tmp[l + num_stars][num_stars + plate_dependent_max_power + 1] = -1
        for p in range(plate_dependent_max_power):
            tmp[l][num_stars + 1 + p] = (stars[l] - plates[0]) ** (p + 1)
            tmp[l + num_stars][num_stars + plate_dependent_max_power + 2 + p] = (stars[l] - plates[1]) ** (p + 1)
        for p in range(plate_independent_max_power - plate_dependent_max_power):
            tmp[l][num_stars + 2 + p + 2 * plate_dependent_max_power] \
                = (stars[l] - plates[0]) ** (p + plate_dependent_max_power + 1)
            tmp[l + num_stars][num_stars + 2 + p + 2 * plate_dependent_max_power] \
                = (stars[l] - plates[1]) ** (p + plate_dependent_max_power + 1)
    return tmp


def make_constraint(num_constraint, num_parameter):
    tmp = np.zeros((num_constraint, num_parameter))
    for i in range(num_constraint):
        tmp[i][i] = 1
    return tmp


p = Plates()
print('insert 3, 5, 1...')
p.insert(-1.0, 1.0, 1)
p.insert(0.0, 2.0, 2)
p.insert(1.0, 3.0, 3)
print('print data')
p.show()
print('delete 3...')
p.delete(2)
print('print data')
p.show()
print('reset id')
p.reset_id()
p.show()


"""
s = Stars()
print('insert 3, 5, 1...')
s.insert(3, 1)
s.insert(5, 2)
s.insert(1, 3)
s.insert(2, 4)
s.insert(9, 5)
print('print data')
s.show()
print('delete 3...')
s.delete(3)
print('print data')
s.show()
print('reset id')
s.reset_id()
s.show()
"""
"""
L = 40

xi = make_stars(L)
xc = make_plate()

for n in range(7):
    for m in range(n + 1):
        # fname = 'a.csv'
        print("n={}, m={}".format(n, m))
        a = make_design_matrix(L, m, n, xi, xc)
        c = make_constraint(8, a.shape[1])
        a = np.concatenate((a, c))
        u, s, v = np.linalg.svd(a)
        # np.savetxt(fname, a, delimiter=',')
        tmp = np.arange(a.shape[1])
        print(s)
        # with open(fname, 'a') as f_handle:
        #    np.savetxt(f_handle, u, delimiter=',')
        #    np.savetxt(f_handle, [tmp], delimiter=',')
        #    np.savetxt(f_handle, v, delimiter=',')
        #    np.savetxt(f_handle, s, delimiter=',')


u, s, v = np.linalg.svd(a)

# print(xi)
# print(xc)
# np.savetxt('a.csv', a, delimiter=',')
"""
