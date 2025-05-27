from fractions import Fraction
import math
import sys

try:
    sys.stdout.reconfigure(encoding='utf-8')
except:
    pass

def invert_matrix(mat):
    """возвращает обратную матрицу"""
    n = len(mat)
    aug = [row[:] + [Fraction(int(i == j)) for j in range(n)] for i, row in enumerate(mat)]
    for i in range(n):
        if aug[i][i] == 0:
            for r in range(i+1, n):
                if aug[r][i] != 0:
                    aug[i], aug[r] = aug[r], aug[i]
                    break
        piv = aug[i][i]
        if piv == 0:
            raise RuntimeError("оббратная матрица не существует")
        for j in range(2*n):
            aug[i][j] /= piv
        for r in range(n):
            if r != i:
                factor = aug[r][i]
                if factor != 0:
                    for j in range(i, 2*n):
                        aug[r][j] -= factor * aug[i][j]
    inv = [row[n:] for row in aug]
    return inv

class ParametricSimplex:
    def __init__(self, A, b_prime, b_double, c):
        # сохраняем данные как дроби
        self.orig_A = [[Fraction(x) for x in row] for row in A]
        self.orig_b_prime = [Fraction(x) for x in b_prime]
        self.orig_b_double = [Fraction(x) for x in b_double]
        self.orig_c = [Fraction(x) for x in c]

        self._setup(self.orig_A, self.orig_b_prime, self.orig_b_double, self.orig_c) #добавление базисных столбцов

    def _setup(self, A, b_prime, b_double, c):
        """подготавливает расширенную матрицу, добавляет свободные переменные и начальный базис."""
        self.A = [row[:] for row in A]
        self.b_prime = b_prime[:]
        self.b_double = b_double[:]
        self.c = c[:]
        self.m = len(self.A)                   # колво ограничений
        self.n = len(self.A[0]) if self.m > 0 else 0  # колва переменных 
        # Добавляем столбцы запасных переменных (единичная матрица m x m)
        for i in range(self.m):
            self.A[i].extend(Fraction(1) if j == i else Fraction(0) for j in range(self.m))
        # Коэффициенты функции цели для запасных переменных = 0
        self.c.extend(Fraction(0) for _ in range(self.m))
        self.total_vars = self.n + self.m
        # Изначальный базис: индексы запасных переменных
        self.basis = list(range(self.n, self.n + self.m))
        self.nonbasis = list(range(self.n))
        # Вычисляем обратную матрицу базисной матрицы B
        self._update_B_inv()

    def _update_B_inv(self):
        """Вычисляет обратную матрицу текущей базисной подматрицы B."""
        B_cols = []
        for bi in self.basis:
            col = [self.A[i][bi] for i in range(self.m)]
            B_cols.append(col)
        B = [list(row) for row in zip(*B_cols)]  # транспонируем список столбцов
        self.B_inv = invert_matrix(B)

    def _simplex(self, b_const):
        """
         симплекс-метод. 
        """
        for _ in range(1000):
            # Вычисляем текущее базисное решение x_B = B_inv * b_const
            xB = [sum(self.B_inv[i][k] * b_const[k] for k in range(self.m)) for i in range(self.m)]
            # Коэффициенты ЦФ для базисных переменных
            cB = [self.c[j] for j in self.basis]
            # Вычисляем вектор y = cB^T * B_inv
            y = [sum(cB[k] * self.B_inv[k][i] for k in range(self.m)) for i in range(self.m)]
            # Выбираем входящую переменную (минимальный индекс с positive reduced cost)
            enter = None
            for j in sorted(self.nonbasis):
                Aj = [self.A[i][j] for i in range(self.m)]
                reduced = self.c[j] - sum(y[k] * Aj[k] for k in range(self.m))
                if reduced > 0:
                    enter = j
                    break
            # Оптимум достигнут, если не нашли входящую переменную
            if enter is None:
                return
            # Находим столбец B_inv * A[:, enter]
            Aj = [self.A[i][enter] for i in range(self.m)]
            col = [sum(self.B_inv[i][k] * Aj[k] for k in range(self.m)) for i in range(self.m)]
            # Правило отношения: находим уходящую переменную
            leave_idx = None
            min_ratio = None
            for i in range(self.m):
                if col[i] > 0:
                    ratio = xB[i] / col[i]
                    if leave_idx is None or ratio < min_ratio or (ratio == min_ratio and self.basis[i] < self.basis[leave_idx]):
                        min_ratio = ratio
                        leave_idx = i
            if leave_idx is None:
                raise RuntimeError("решение неограничено")
            # Обмен базисов: входящая -> уходящая
            leaving_var = self.basis[leave_idx]
            self.basis[leave_idx] = enter
            self.nonbasis.remove(enter)
            self.nonbasis.append(leaving_var)
            # Обновляем обратную матрицу
            self._update_B_inv()
        raise RuntimeError("более 1000итераций")

    def _compute_interval(self):
        """
        Вычисляет допустимый интервал λ (L, U) для текущего базиса.
        Возвращает (L, U, alpha, beta) или None, если базис недопустим.
        """
        # alpha = B_inv * b_prime, beta = B_inv * b_double
        alpha = [sum(self.B_inv[i][k] * self.b_prime[k] for k in range(self.m)) for i in range(self.m)]
        beta  = [sum(self.B_inv[i][k] * self.b_double[k] for k in range(self.m)) for i in range(self.m)]
        L, U = -math.inf, math.inf
        for i in range(self.m):
            if beta[i] > 0:
                val = -alpha[i] / beta[i]
                if val > L:
                    L = val
            elif beta[i] < 0:
                val = -alpha[i] / beta[i]
                if val < U:
                    U = val
            else:
                # Если beta[i] == 0, то alpha[i] должно быть >= 0
                if alpha[i] < 0:
                    return None  # Базис недопустим
        # Проверяем корректность интервала (L ≤ U)
        if L > U:
            return None
        return (L, U, alpha, beta)

    def _dual_enter(self, leave_idx):
        """Для заданной уходящей строки выбирает входящую переменную."""
        # Строка (B_inv * A)[leave_idx]
        row = [sum(self.B_inv[leave_idx][k] * self.A[k][j] for k in range(self.m)) for j in range(self.total_vars)]
        # Кандидаты: небазисные с отрицательным элементом
        candidates = [j for j in self.nonbasis if row[j] < 0]
        return None if not candidates else min(candidates)

    def _find_initial_lambda(self):
        """
        Находит допустимое значение λ для начального базиса.
        Если не существует, выбрасывает понятное исключение.
        
        """
        lam_min, lam_max = -math.inf, math.inf
        for i in range(self.m):
            alpha_i = self.orig_b_prime[i]
            beta_i  = self.orig_b_double[i]
            if beta_i > 0:
                lam_min = max(lam_min, -alpha_i / beta_i)
            elif beta_i < 0:
                lam_max = min(lam_max, -alpha_i / beta_i)
            else:
                if alpha_i < 0:
                    raise RuntimeError("Начальный базис не имеет допустимого λ (противоречие)")
        if lam_min > lam_max:
            raise RuntimeError("Невозможно найти допустимое λ для начального базиса")
        # Если λ может быть любым, выбираем 0
        if lam_min == -math.inf and lam_max == math.inf:
            return Fraction(0)
        if lam_min == -math.inf:
            return lam_max
        if lam_max == math.inf:
            return lam_min
        return Fraction(0) if lam_min <= 0 <= lam_max else lam_min

    def solve(self):
        """
        Главный метод: решает параметрическую задачу.
        Возвращает список кортежей (L, U, alpha, beta, basis) для каждого интервала λ.
        """
        # Находим начальное λ и решаем симплекс-метод для него
        lam0 = self._find_initial_lambda()
        b0 = [self.b_prime[i] + lam0 * self.b_double[i] for i in range(self.m)]
        self._simplex(b0)

        intervals = []
        # Первый интервал для исходного базиса
        iv = self._compute_interval()
        if iv is None:
            return intervals
        L0, U0, alpha0, beta0 = iv
        intervals.append((L0, U0, alpha0, beta0, self.basis.copy()))

        # Шаги по λ вправо (увеличение)
        currL, currU, currAlpha, currBeta = L0, U0, alpha0, beta0
        while currU != math.inf:
            # Находим переменную, уходящую при λ = currU
            leave_idx = None
            for i in range(self.m):
                if currBeta[i] < 0 and -currAlpha[i] / currBeta[i] == currU:
                    leave_idx = i
                    break
            if leave_idx is None:
                break
            # Находим входящую по двойственному правилу
            enter = self._dual_enter(leave_idx)
            if enter is None:
                break
            # Обмен базисных переменных
            leaving_var = self.basis[leave_idx]
            self.basis[leave_idx] = enter
            self.nonbasis.remove(enter)
            self.nonbasis.append(leaving_var)
            self._update_B_inv()
            # Новый интервал
            iv2 = self._compute_interval()
            if iv2 is None:
                break
            L2, U2, alpha2, beta2 = iv2
            intervals.append((currU, U2, alpha2, beta2, self.basis.copy()))
            currL, currU, currAlpha, currBeta = L2, U2, alpha2, beta2

        # Шаги по λ влево (уменьшение)
        # Восстанавливаем исходный базис
        self._setup(self.orig_A, self.orig_b_prime, self.orig_b_double, self.orig_c)
        lam0 = self._find_initial_lambda()
        b0 = [self.b_prime[i] + lam0 * self.b_double[i] for i in range(self.m)]
        self._simplex(b0)
        iv = self._compute_interval()
        if iv is None:
            return intervals
        L0, U0, alpha0, beta0 = iv
        currL, currU, currAlpha, currBeta = L0, U0, alpha0, beta0
        while currL != -math.inf:
            leave_idx = None
            for i in range(self.m):
                if currBeta[i] > 0 and -currAlpha[i] / currBeta[i] == currL:
                    leave_idx = i
                    break
            if leave_idx is None:
                break
            enter = self._dual_enter(leave_idx)
            if enter is None:
                break
            leaving_var = self.basis[leave_idx]
            self.basis[leave_idx] = enter
            self.nonbasis.remove(enter)
            self.nonbasis.append(leaving_var)
            self._update_B_inv()
            iv2 = self._compute_interval()
            if iv2 is None:
                break
            L2, U2, alpha2, beta2 = iv2
            intervals.insert(0, (L2, currL, alpha2, beta2, self.basis.copy()))
            currL, currU, currAlpha, currBeta = L2, U2, alpha2, beta2

        return intervals

# пример
if __name__ == '__main__':

    c = [3, 0, 0, -3, 2] 
    A = [
        [0, 1, 0, 3, -10],   
        [1, 0, 0, -1, 13],
        [0, 0, 1, 7, 6]    
    ]
    b_prime = [27, 117, 151]     
    b_double = [-18, 18, -18]  

    ps = ParametricSimplex(A, b_prime, b_double, c)
    intervals = ps.solve()

    if intervals:
        # Проверка: нет решений до первого интервала
        if intervals[0][0] != -math.inf:
            print(f"[-\u221E, {intervals[0][0]}): решений нет\n")

        for (L, U, alpha, beta, basis) in intervals:
            # Вычисляем коэффициенты функции цели: F = c0 + c1*λ
            c0 = Fraction(0)
            c1 = Fraction(0)
            for i, var in enumerate(basis):
                if var < len(ps.orig_c):
                    c0 += ps.orig_c[var] * alpha[i]
                    c1 += ps.orig_c[var] * beta[i]
            # Формируем решение X
            X_expr = []
            print(basis)
            for j in range(ps.total_vars):
                if j in basis:
                    idx = basis.index(j)
                    a_j = alpha[idx]
                    b_j = beta[idx]
                else:
                    a_j = Fraction(0)
                    b_j = Fraction(0)
                if b_j >= 0:
                    X_expr.append(f"{a_j} + {b_j}*λ" if b_j != 0 else f"{a_j}")
                else:
                    X_expr.append(f"{a_j} - {-b_j}*λ")

            # Строковые представления границ
            L_str = "-\u221E" if L == -math.inf else str(L)
            U_str = "\u221E" if U == math.inf else str(U)

            print(f"[{L_str}, {U_str}]")
            print(f"F = {c0}" + (f" + {c1}*λ" if c1 >= 0 else f" - {-c1}*λ"))
            print("X = [" + ", ".join(X_expr) + "]\n")

        # Проверка: нет решений после последнего интервала
        if intervals[-1][1] != math.inf:
            print(f"({intervals[-1][1]}, \u221E]: решений нет\n")
