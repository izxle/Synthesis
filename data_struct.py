from periodictable import elements
from common_names import parse_formula
import copy

parse_element = elements.symbol


class Atom:
    def __init__(self, element, mass, moles, wt_percentage, at_percentage, precursor=None):
        self.element = element
        self.mass = mass
        self.moles = moles
        self.wt_percentage = wt_percentage
        self.at_percentage = at_percentage
        self.precursor = precursor

    def set_precursor(self, precursor):
        assert self.element in precursor.atoms, f"{precursor} is not a valid precursor for {self.element}"
        self.precursor = precursor

    def __format__(self, format_spec):
        return str(self)

    def __str__(self):
        return str(self.element)


class Compound:
    def __init__(self, formula, mass, moles, purity):
        self.formula = formula
        self.atoms = formula.atoms
        self.mass = mass
        self.moles = moles
        self.purity = purity

    def copy(self):
        return Compound(self.formula, self.mass, self.moles, self.purity)

    def __format__(self, format_spec):
        return str(self)

    def __str__(self):
        return str(self.formula)

    def __eq__(self, other):
        return self.atoms == other.atoms


class Synthesis:
    def __init__(self, composition: dict, mass: float, precursors: dict,
                 percentage: str='atomic', experimental=None):
        assert percentage in ['atomic', 'weight'], f"percentage must be 'atomic' or 'weight', not {percentage}"
        self.composition = dict()
        self.atoms = dict()
        self.experimental = experimental
        self.mass = mass
        self.percentage = percentage
        self.set_atoms(composition)
        self.set_precursors(precursors)
        if experimental:
            self._calc_experimental()

    def set_precursors(self, precursors):
        for prec, *info in precursors.items():
            self.add_precursor(prec, info)
        self.check_precursors()

    def _parse_percentages(self, composition):
        atoms = {parse_element(e): {self.percentage: p / 100} for e, p in composition.items()}

        if self.percentage == 'atomic':
            total_mass = sum(element.mass * percentage['atomic']
                             for element, percentage in atoms.items())
            for element, percentage in atoms.items():
                atoms[element]['weight'] = element.mass * percentage['atomic'] / total_mass

        elif self.percentage == 'weight':
            total_moles = sum(percentage['weight'] / element.mass
                              for element, percentage in atoms.items())
            for element, percentage in atoms.items():
                atoms[element]['atomic'] = (percentage['weight'] / element.mass) / total_moles

        return atoms

    def set_atoms(self, composition: dict):
        atoms = self._parse_percentages(composition)
        self.moles = sum(self.mass * percentage['weight'] / element.mass
                         for element, percentage in atoms.items())

        for element, percentage in atoms.items():
            self.composition[element] = percentage
            mass = self.mass * percentage['weight']
            moles = self.moles * percentage['atomic']
            self.atoms[element] = Atom(element=element,
                                       mass=mass,
                                       moles=moles,
                                       at_percentage=percentage['atomic'],
                                       wt_percentage=percentage['weight'])

    def parse_info(self, info: list):
        res = dict()
        for val in info:
            if '%' in val:
                res['purity'] = float(val.strip()[:-1]) / 100
        return res

    def add_precursor(self, formula, info):
        info = self.parse_info(info)
        purity = info.get('purity', 1)

        precursor = parse_formula(formula)
        element = next(elm for elm in self.composition if elm in precursor.atoms)
        n = precursor.atoms[element]
        moles = self.atoms[element].moles / n
        mass = moles * precursor.mass / purity
        precursor = Compound(formula=precursor,
                             purity=purity,
                             mass=mass,
                             moles=moles)
        self.atoms[element].set_precursor(precursor)

    def check_precursors(self):
        # restricted to only one precursor per metal
        assert all(a.precursor for a in self.atoms.values()), f"found no precursors for element '{a.element}'"
        # for element in self.composition:
        #     prec = [prec for prec in self.precursors if element in prec.atoms]
        #     assert len(prec) == 1, f'There must be 1 precursor for element {element}, found {len(prec)}: {prec}'

    def _calc_experimental(self):
        # create copy of self.atoms
        atoms = copy.deepcopy(self.atoms)
        total_mass = 0
        total_moles = 0
        for formula, mass in self.experimental.items():
            precursor = parse_formula(formula)
            element = next(elm for elm in self.composition if elm in precursor.atoms)
            prev_precursor = self.atoms[element].precursor
            assert prev_precursor == precursor, f"precursor for {element} does not match"
            purity = prev_precursor.purity
            moles = mass * purity / precursor.mass
            new_precursor = Compound(formula=precursor,
                                     purity=purity,
                                     mass=mass,
                                     moles=moles)
            atoms[element].set_precursor(new_precursor)
            total_mass += mass * purity * element.mass / precursor.mass
            total_moles += moles

        self.experimental['mass'] = total_mass
        self.experimental['moles'] = total_moles

        for element, atom in atoms.items():
            precursor = atom.precursor
            n = precursor.atoms[element]
            moles = precursor.moles * n
            mass = moles * element.mass

            atoms[element].mass = mass
            atoms[element].moles = moles
            atoms[element].at_percentage = moles / total_moles
            atoms[element].wt_percentage = mass / total_mass
        self.experimental['atoms'] = atoms


    def __str__(self):
        res = 'Theoretical\n'
        res += f'Core:      {self.mass / 1e3:.1f} mg\t{self.moles / 1e3:.1f} mmols\n'
        res += 'mass [mg]\tmoles [mmol]\n'
        res += 'Element   at%   wt%   mass0   massf     moles   precursor       mass\n'
        # TODO: add reation
        atoms = (f'   {a:3}   {a.at_percentage * 100:5.1f} {a.wt_percentage * 100:5.1f}  ' +
                 f'{a.mass * 1e3:6.2f}  {a.mass * 1e3:6.2f}  ' +
                 f'{a.moles * 1e3:8.5f}   {a.precursor:20}   {a.precursor.mass * 1e3:6.3f}'
                 for a in self.atoms.values())
        res += '\n'.join(atoms) + '\n'
        
        if self.experimental:
            res += f'Experimental\nCore:      {self.experimental["mass"] / 1e3:.1f} mg\t'
            res += f'{self.experimental["moles"] / 1e3:.1f} mmols\n'
            res += 'mass [mg]\tmoles [mmol]\n'
            res += 'Element   at%   wt%   mass0   massf     moles   precursor       mass\n'
            # TODO: add reation
            atoms = (f'   {a:3}   {a.at_percentage * 100:5.1f} {a.wt_percentage * 100:5.1f}  ' +
                     f'{a.mass * 1e3:6.2f}  {a.mass * 1e3:6.2f}  ' +
                     f'{a.moles * 1e3:8.5f}   {a.precursor:20}   {a.precursor.mass * 1e3:6.3f}'
                     for a in self.experimental["atoms"].values())
            res += '\n'.join(atoms)

        return res
