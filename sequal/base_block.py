
class BaseBlock:
    def __init__(self, value, position, branch=False, mass=None):
        self.value = value
        self.position = position
        self.branch = branch
        self.mass = mass
        self.extra = None

    def __str__(self):
        return self.value

    def __repr__(self):
        return self.value