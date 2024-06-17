from collections import deque
from termcolor import colored

"""
It is known that Clifford QCA on a line fall into one of 3 classes:
1. Glider
2. Fractal
3. Periodic
This is a simple tool to simulate QCA dynamics on a line.
"""

I, X, Y, Z = 0, 1, 2, 3
# Pauli Products; for now, not tracking phases
PP_dict = {
    (I, I): I, (I, X): X, (I, Y): Y, (I, Z): Z,
    (X, I): X, (X, X): I, (X, Y): Z, (X, Z): Y,
    (Y, I): Y, (Y, X): Z, (Y, Y): I, (Y, Z): X,
    (Z, I): Z, (Z, X): Y, (Z, Y): X, (Z, Z): I,
}

FILLED_CIRCLE_CHR = '\u25CF'
FILLED_BOX_CHR = '\u2588' # '\u25A0', '\u25A3', '\u25A6' also good
X_MARK_CHR = '\u2716' # '\u2573' also good

DARK_COLOR_SCHEME = ['black', 'light_red', 'light_blue', 'green']
LIGHT_COLOR_SCHEME = ['white', 'yellow', 'blue', 'green']

COMMON_PRINT_SYMBOL = FILLED_BOX_CHR
COLOR_SCHEME = LIGHT_COLOR_SCHEME

# 1D translationally invariant Clifford QCA
class CliffordQCA:
    print_symbols_v0 = {I: ' ', X: 'x', Y: 'y', Z: 'z'}
    print_symbols_v1 = {
        I: colored(' ', 'white'),
        X: colored('x', 'red'),
        Y: colored('y', 'yellow'),
        Z: colored('z', 'green'),
    }
    print_symbols_v2 = {
        I: colored(COMMON_PRINT_SYMBOL, COLOR_SCHEME[0]),
        X: colored(COMMON_PRINT_SYMBOL, COLOR_SCHEME[1]),
        Y: colored(COMMON_PRINT_SYMBOL, COLOR_SCHEME[2]),
        Z: colored(COMMON_PRINT_SYMBOL, COLOR_SCHEME[3]),
    }
    # static methods
    def derive_y_rule(x_rule, z_rule):
        left_op = PP_dict[(x_rule[0], z_rule[0])]
        center_op = PP_dict[(x_rule[1], z_rule[1])]
        right_op = PP_dict[(x_rule[2], z_rule[2])]
        return (left_op, center_op, right_op)

    def __init__(self, l, n, x_rule, z_rule, state=None, debug=False):
        self.l = l # TODO: assert l < (n - 1) // 2 (for robust indexing)
        self.n = n
        self.middle_qubit_idx = (n - 1) // 2 # NB: in the even case, pick left
        self.x_rule = x_rule
        self.z_rule = z_rule
        self.state = state if state else [I] * n
        self.curr_time_step = 0
        self.history = []
        self.debug = debug
        self.all_ca_rules = self.derive_all_ca_rules() # TODO: move this/refactor

    # derive rule for qubit in center of lattice
    def derive_op_update_rule(self, op):
        rule = (I, I, I)
        if op == I:
            pass
        if op == X:
            rule = self.x_rule
        elif op == Y:
            rule = CliffordQCA.derive_y_rule(self.x_rule, self.z_rule)
        elif op == Z:
            rule = self.z_rule
        else:
            RuntimeError('Invalid op: {}'.format(op))

        return rule

    def convert_rule_to_pauli(self, rule, offset=0):
        res = deque([I] * self.n)

        center_idx = self.middle_qubit_idx
        left_idx = center_idx - self.l
        right_idx = center_idx + self.l

        res[left_idx], res[center_idx], res[right_idx] = rule
        res.rotate(offset)
        return res

    # TODO: better name for this...
    def combine_evo(self, left_rule, center_rule, right_rule):
        left_evo = self.convert_rule_to_pauli(left_rule, offset=-self.l)
        center_evo = self.convert_rule_to_pauli(center_rule)
        right_evo = self.convert_rule_to_pauli(right_rule, offset=self.l)

        # int_res = left_evo * center_evo
        int_res = [PP_dict[(left_op, right_op)]
                   for left_op, right_op in zip(left_evo, center_evo)]
        # combined_res = left_evo * center_evo * right_evo
        return [PP_dict[(left_op, right_op)]
                   for left_op, right_op in zip(int_res, right_evo)]

    # derive rules for qubit in center of lattice; shift result tuple for other sites
    def derive_all_ca_rules(self):
        possibilites = [(op1, op2, op3)
                        for op1 in [I, X, Y, Z]
                        for op2 in [I, X, Y, Z]
                        for op3 in [I, X, Y, Z]]

        res = {}
        for p in possibilites:
            left_op, center_op, right_op = p

            left_op_rule = self.derive_op_update_rule(left_op)
            center_op_rule = self.derive_op_update_rule(center_op)
            right_op_rule = self.derive_op_update_rule(right_op)

            # updated_center_op = left_op_rule[2] * center_op_rule[1] * right_op_rule[0]
            l_part, c_part, r_part = left_op_rule[2], center_op_rule[1], right_op_rule[0]
            l_part_times_c_part = PP_dict[(l_part, c_part)]
            res[p] = PP_dict[(l_part_times_c_part, r_part)]

            if self.debug:
                evolved_pauli = self.combine_evo(left_op_rule, center_op_rule, right_op_rule)
                print(p, evolved_pauli, res[p])

        return res

    # TODO: choose a better name
    def partition_lattice(self):
        all_tups = []
        for center_idx in range(self.n):
            left_idx = center_idx - self.l
            right_idx = (center_idx + self.l) % self.n
            all_tups.append(
                (self.state[left_idx], self.state[center_idx], self.state[right_idx]))
        return all_tups

    def evolve_ca(self, time_steps=1):
        for _ in range(time_steps):
            self.history.append(self.state)
            self.state = [self.all_ca_rules[tup] for tup in self.partition_lattice()]
        self.curr_time_step += time_steps

    def __str__(self):
        print_symbols = CliffordQCA.print_symbols_v2
        buffer = ''
        for state in self.history + [self.state]:
            state_str = ''.join([print_symbols[symbol] for symbol in state])
            buffer = state_str + '\n' + buffer
        return buffer

#####################################################
# Glider Examples
#####################################################
# n = 31
# starting_state = [I]*14 + [X, Z] + [I]*15 # glider
# test = CliffordQCA(1, n, (I, Z, I), (Z, X, Z), starting_state)

# time_steps = 62
# test.evolve_ca(time_steps)
# print(test)

#####################################################
# Fractal Examples
#####################################################
# NB: Interestingly, if time_steps <= side_sz, the dyanmics look clean, BUT
# if time_steps > side_sz, the dynamics look unstructured. It sort of seems
# like there's some reflections muddying everything i.e. feels like wave-like
# dynamics/behavior. Could be totally off, but curious about this.

side_sz = 50 # max size that we can visualize in terminal is ~325
n = 2 * side_sz + 1
starting_state = [I] * side_sz + [X] + [I] * side_sz
fractal_qca = CliffordQCA(1, n, (X, Y, X), (I, X, I), starting_state)

time_steps = 50
fractal_qca.evolve_ca(time_steps)
print(fractal_qca)

#####################################################
# Periodic Examples
#####################################################
# n = 31
# starting_state = [I]*15 + [X] + [I]*15
# test = CliffordQCA(1, n, (I, Z, I), (Z, X, Z), starting_state)

# time_steps = 123
# test.evolve_ca(time_steps)
# print(test)
