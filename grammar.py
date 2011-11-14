class Grammar:
    """ Grammar is just a dictionary that represents
        a grammar in Backus-Naur Form and that uses
        module when acceding a rule.
    """
    def __init__(self, bnf):
        assert isinstance(bnf, dict)
        self._bnf = bnf

    def __getitem__(self, args):
        if isinstance(args, str):
            # Se llamo con un solo argumento
            return self._bnf[args]
        else: 
            assert len(args) == 2
            rule, option = args
            size_rule = len(self._bnf[rule])
            return self._bnf[rule][option % size_rule]
