class TaxClade:
    def __init__(self, name):
        self.children = {}
        self.name, self.father = name, None
        self.n_children = 0

    def update_child_count(self):
        self.n_children += 1
        cl = self.father
        while cl:
            cl.n_children += 1
            cl = cl.father

    def add_child(self, name):
        new_clade = TaxClade(name)
        self.children[name] = new_clade
        self.update_child_count()
        new_clade.father = self
        return new_clade

    def get_terminals(self):
        terms = []
        if not self.children:
            return [self]
        for c in list(self.children.values()):
            terms += c.get_terminals()
        return terms

    def get_full_name(self):
        fullname = [self.name]
        cl = self.father
        while cl:
            fullname = [cl.name] + fullname
            cl = cl.father
        return "|".join(fullname[1:])


class TaxTree:
    def __init__(self, mpa_pkl):
        self.root = TaxClade("root")
        self.all_clades = {}
        clades_txt = ((l.strip().split('|')[:-1])
                      for l, n in list(mpa_pkl['taxonomy'].items()))
        for clade in clades_txt:
            father = self.root
            for clade_lev in clade:
                if not clade_lev in father.children:
                    father.add_child(clade_lev)
                    self.all_clades[clade_lev] = father.children[clade_lev]
                father = father.children[clade_lev]

    def get_lineage(self, lineage):
        cl = self.root
        for clade in lineage.split('|'):
            cl = cl.children[clade]
        return cl
