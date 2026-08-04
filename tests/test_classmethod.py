class A:
    @classmethod
    def __call__(cls, x):
        return cls.helper(x)

    @classmethod
    def helper(cls, x):
        return x * 2


print(A()(5))
