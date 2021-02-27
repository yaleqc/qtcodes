from abc import abstractmethod, ABCMeta


class ExMixin:
    def dostuff(self) -> None:
        print("wowow")


class Base:
    def __init__(self, a):
        self.a = 2 * a


class MainStuff(Base):
    def __init__(self, a, b):
        super().__init__(a)
        self.b = b


class MixedClass(ExMixin, MainStuff):
    pass


hm = MixedClass(2, 3)

print(hm.a)
print(hm.b)
hm.dostuff()
