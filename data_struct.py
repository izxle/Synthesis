from periodictable import elements
from common_names import parse_formula
from copy import deepcopy
from typing import Tuple
from math import pi

parse_element = elements.symbol


class SimpleReagent:
    def __init__(self, name, mass=None, moles=None, volume=None):
        self.name = name
        self.mass = mass
        self.moles = moles
        self.volume = volume

    def __format__(self, format_spec):
        return f'{str(self):{format_spec}}'

    @property
    def quantity(self):
        if self.mass is not None:
            res = f'{self.mass:.1f} mg'
        elif self.volume is not None:
            res = f'{self.volume:.3f} mL'
        else:
            res = str(self)
        return res

    def __str__(self):
        return str(self.name)

    def __eq__(self, other):
        return self.name == other.name


class Atom(SimpleReagent):
    def __init__(self, element, mass, moles, wt_percentage, at_percentage, precursor=None):
        super().__init__(name=element, mass=mass, moles=moles)
        self.element = element
        self.wt_percentage = wt_percentage
        self.at_percentage = at_percentage
        self.precursor = precursor

    def set_precursor(self, precursor):
        assert self.element in precursor.atoms, f"{precursor} is not a valid precursor for {self.element}"
        self.precursor = precursor


class Compound(SimpleReagent):
    def __init__(self, formula, mass, moles, purity=1):
        super().__init__(name=formula, mass=mass, moles=moles)
        self.formula = formula
        self.atoms = formula.atoms
        self.purity = purity

    def copy(self):
        return Compound(self.name, self.mass, self.moles, self.purity)

    def __eq__(self, other):
        return self.atoms == other.atoms


class Experimental:
    def __init__(self, atoms, precursors: dict, support: Tuple[str, float]=None):
        """

        :param atoms: every Atom must have a precursor
        :param precursors:
        :param support:
        """
        self.mass = 0
        self.catalyst_mass = 0
        self.catalyst_moles = 0
        self.wt_percentage = 1
        self.atoms = atoms
        self.set_precursors(precursors)
        self.set_support(*support)

    def set_support(self, name, mass):
        support = SimpleReagent(name=name, mass=mass)
        self.support = support
        self.mass += mass
        self.wt_percentage = self.catalyst_mass / self.mass

    def set_precursors(self, precursors: dict):
        atoms = self.atoms

        # TODO: check in readers
        assert all(parse_element(e) in atoms for e in precursors), 'All elements must be specified in precursors.'

        total_moles = 0
        total_mass = 0
        for element_symbol, mass in precursors.items():
            element = parse_element(element_symbol)
            precursor_formula = atoms[element].precursor.formula

            purity = atoms[element].precursor.purity
            moles = mass * purity / precursor_formula.mass

            precursor = Compound(formula=precursor_formula,
                                 mass=mass,
                                 moles=moles,
                                 purity=purity)

            atoms[element].set_precursor(precursor)

            element_moles = moles * precursor_formula.atoms[element]
            element_mass = element_moles * element.mass

            total_moles += element_moles
            total_mass += element_mass

            atoms[element].mass = element_mass
            atoms[element].moles = element_moles

        self.catalyst_mass = total_mass
        self.catalyst_moles = total_moles
        self.mass = total_mass

        # calculating element weight and atomic percentages in final catalyst
        for element, atom in atoms.items():
            atoms[element].at_percentage = atom.moles / total_moles
            atoms[element].wt_percentage = atom.mass / total_mass


class ThinFilm:
    def __init__(self, load: float, diameter: float,
                 support_wt_percentage: float=0, PGM_wt_percentage: float=1.0):
        self.load = load
        self.support_wt_percentage = support_wt_percentage
        self.PGM_wt_percentage = PGM_wt_percentage
        self.set_area(diameter=diameter)

        PGM_mass = load * self.area
        catalyst_mass = PGM_mass / PGM_wt_percentage
        total_mass = catalyst_mass * (2 - support_wt_percentage)
        self.PGM_mass = PGM_mass
        self.total_mass = total_mass

    def ink_from_sample(self, sample_volume, ink_volume=None, catalyst_mass=None):
        ink_concentration = self.total_mass / sample_volume
        if ink_volume:
            catalyst_mass = ink_concentration * ink_volume
        elif catalyst_mass:
            ink_volume = catalyst_mass / ink_concentration
        else:
            raise ValueError('must provide `ink_volume` or `catalyst_mass`.')
        ink = Ink(mass=catalyst_mass, volume=ink_volume, concentration=ink_concentration)
        self.ink = ink
        return ink

    def set_area(self, area=None, diameter=None):
        if area:
            self.area = area
        elif diameter:
            # diameter in cm
            self.diameter = diameter
            self.area = pi * (self.diameter / 2)**2
        else:
            raise ValueError('Must provide area or diameter')

    def __str__(self):
        res = 'Thin Film\n' \
              'Load     Area    mass    PGM_mass\n' \
              'ug/cm^2  cm^2    ug      ug\n' \
              f'{self.load:7.1f}  {self.area:6.3f}  {self.total_mass:6.3f}  {self.PGM_mass:7.3f}\n  '
        return res


class Ink:
    def __init__(self, mass, volume, concentration):
        self.mass = mass
        self.volume = volume
        self.concentration = concentration

    def __str__(self):
        res = 'Ink\n' \
              'mass    volume  Conc\n' \
              'mg      mL      mg/mL\n' \
              f'{self.mass:6.4f}  {self.volume:6.1f}  {self.concentration:5.3f}\n'
        return res


class Synthesis:
    def __init__(self, composition: dict, mass: float, precursors: dict, percentage: str='atomic',
                 support: tuple=None,
                 other_reagents: dict = None,
                 experimental: dict=None,
                 ink: dict=None):
        assert percentage in ['atomic', 'weight'], f"percentage must be 'atomic' or 'weight', not {percentage}"
        self.composition = dict()
        self.atoms = dict()
        self.experimental = dict()
        self.ink = dict()
        self.mass = mass
        self.percentage = percentage
        self.set_atoms(composition)
        self.set_precursors(precursors)
        # non-essential information
        self.set_support(support)
        self.set_other_reagents(other_reagents)
        if experimental:
            self.set_experimental(experimental)
            if ink:
                self.set_ink(ink)

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
        if isinstance(info, str):
            info = [info]
        for val in info:
            if '%' in val:
                res['purity'] = float(val.strip()[:-1]) / 100
            elif ':' in val:
                this, ref = val.split(':')
                res['this'] = dict(value=float(this[:-1]),
                                   unit=this[-1])
                res['ref'] = dict(value=float(ref[:-1]),
                                  unit=ref[-1])
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

    def set_experimental(self, experimental: dict):
        # create copy of self.atoms
        atoms = deepcopy(self.atoms)
        self.experimental = Experimental(atoms=atoms, **experimental)

    def set_ink(self, ink: dict):
        PGM = (parse_element('Pt'), parse_element('Pd'))
        PGM_wt_percentage = sum(a.wt_percentage
                                for a in self.experimental.atoms.values()
                                if a.element in PGM)
        support_wt_percentage = self.experimental.support.mass / self.experimental.mass

        load = ink['load']
        diameter = ink['diameter']
        film = ThinFilm(load=load, diameter=diameter,
                        support_wt_percentage=support_wt_percentage,
                        PGM_wt_percentage=PGM_wt_percentage)
        catalyst_mass = ink.get('catalyst_mass')
        ink_volume = ink.get('ink_volume')
        sample_volume = ink['sample_volume']
        ink = film.ink_from_sample(catalyst_mass=catalyst_mass, ink_volume=ink_volume, sample_volume=sample_volume)

        self.film = film
        self.ink = ink

    def set_support(self, support):
        name, percentage = support
        mass = self.mass * percentage / (100 - percentage)
        self.support = SimpleReagent(name, mass)

    def set_other_reagents(self, other_reagents: dict):
        reagents = []
        for name, ratios in other_reagents.items():
            info = self.parse_info(ratios)
            if info['ref']['unit'] == 'w':
                ref = self.mass
            elif info['ref']['unit'] == 'a':
                ref = self.moles
            else:
                raise NotImplementedError("Only 'a' and 'w' are implemented")

            ref /= info['ref']['value']

            new_val = info['this']['value'] * ref

            mass = None
            moles = None
            volume = None
            if info['this']['unit'] == 'v':
                volume = new_val
            elif info['this']['unit'] == 'w':
                mass = new_val
            else:
                raise NotImplementedError("Only 'v' and 'w' are implemented")
            r = SimpleReagent(name=name, mass=mass, volume=volume, moles=moles)
            reagents.append(r)
        self.other_reagents = reagents

    def __str__(self):
        res = 'Theoretical\n'
        res += f'Core:      {self.mass:.1f} mg\t{self.moles:.1f} mmols\n'
        res += 'mass [mg]\tmoles [mmol]\n'
        res += 'Element    at%   wt%   mass0   massf     moles   precursor         mass\n'
        # TODO: add reaction
        atoms = (f'   {a:3}   {a.at_percentage * 100:5.1f} {a.wt_percentage * 100:5.1f}  ' +
                 f'{a.mass:6.2f}  {a.mass:6.2f}  ' +
                 f'{a.moles:8.5f}   {a.precursor:15}   {a.precursor.mass:6.3f}'
                 for a in self.atoms.values())
        res += '\n'.join(atoms) + '\n'
        res += 'Other reagents\n'
        reagents = (f'{r.name:15} {r.quantity}'
                    for r in self.other_reagents)
        res += '\n'.join(reagents) + '\n'
        if self.support:
            res += f'{self.support:16} {self.support.quantity}\n'
        
        if self.experimental:
            # TODO: set __str__ to Experimental
            res += f'\nExperimental\n' \
                   f'Total:  {self.experimental.mass:.1f} mg    ' \
                   f'Core:  {self.experimental.catalyst_mass:.1f} mg\t' \
                   f'{self.experimental.catalyst_moles:.4f} mmols  ' \
                   f'{self.experimental.wt_percentage * 100:.1f} %wt\n'
            res += 'mass [mg]\tmoles [mmol]\n'
            res += 'Element   at%   wt%   mass0   massf     moles   precursor       mass\n'
            atoms = (f'   {a:3}   {a.at_percentage * 100:5.1f} {a.wt_percentage * 100:5.1f}  ' +
                     f'{a.mass:6.2f}  {a.mass:6.2f}  ' +
                     f'{a.moles:8.5f}   {a.precursor:15}   {a.precursor.mass:6.3f}'
                     for a in self.experimental.atoms.values())
            res += '\n'.join(atoms) + '\n'

            if self.film:
                res += str(self.film)
                res += str(self.ink)


        return res
