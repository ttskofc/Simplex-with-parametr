from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from v5 import ParametricSimplex
from fractions import Fraction
import math

app = FastAPI()

# Настройка CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

class ProblemData(BaseModel):
    matrix: list[list[float]]    # A (ограничения)
    cost: list[float]           # c (целевая функция)
    b: list[float]              # b'
    b_double: list[float]       # b'' (b_double)

def convert_to_fractions(data):
    """Конвертирует float в Fraction для точных вычислений."""
    return [
        [Fraction(str(x)) for x in row]
        if isinstance(row, list)
        else Fraction(str(row))
        for row in data
    ]

def format_fraction(value):
    """Форматирует Fraction в строку с дробью или целым числом"""
    if value == -math.inf:
        return "-∞"
    if value == math.inf:
        return "∞"
    if isinstance(value, Fraction):
        if value.denominator == 1:
            return f"{value.numerator}"
        return f"{value.numerator}/{value.denominator}"
    return str(value)

@app.post("/solve-assignment")
async def solve_problem(data: ProblemData):
    try:
        # Конвертация данных
        A = convert_to_fractions(data.matrix)
        c = convert_to_fractions(data.cost)
        b_prime = convert_to_fractions(data.b)
        b_double = convert_to_fractions(data.b_double)

        # Решение задачи
        ps = ParametricSimplex(A, b_prime, b_double, c)
        intervals = ps.solve()

        lines = []
        # Если до первого интервала нет решений
        output_blocks = []

        # До первого интервала
        if intervals and format_fraction(intervals[0][0]) != "-∞":
            output_blocks.append(f"[-∞, {format_fraction(intervals[0][0])})\nрешений нет")

        for interval in intervals:
            L, U, alpha, beta, basis = interval
            block_lines = []

            # Интервал
            block_lines.append(f"[{format_fraction(L)}, {format_fraction(U)}]")

            # Целевая функция
            c0 = sum(c[var] * alpha[i] for i, var in enumerate(basis) if var < len(c))
            c1 = sum(c[var] * beta[i] for i, var in enumerate(basis) if var < len(c))
            c0_str = format_fraction(c0)
            c1_str = format_fraction(abs(c1))
            block_lines.append(f"F = {c0_str}" + (f" + {c1_str}λ" if c1 > 0 else f" - {c1_str}λ" if c1 < 0 else ""))

            # X
            X_expr = []
            for j in range(ps.total_vars):
                if j in basis:
                    idx = basis.index(j)
                    a_j, b_j = alpha[idx], beta[idx]
                else:
                    a_j, b_j = Fraction(0), Fraction(0)

                if a_j == 0 and b_j == 0:
                    expr = "0"
                else:
                    a_str = format_fraction(a_j)
                    b_str = format_fraction(abs(b_j))
                    if b_j == 0:
                        expr = a_str
                    elif b_j > 0:
                        expr = f"{a_str} + {b_str}λ"
                    else:
                        expr = f"{a_str} - {b_str}λ"
                X_expr.append(expr)

            # Удалить до 3 последних "0"
            count = 0
            for i in range(len(X_expr) - 1, -1, -1):
                if X_expr[i] == "0":
                    X_expr.pop(i)
                    count += 1
                    if count == 3:
                        break

            block_lines.append(f"X = [{', '.join(X_expr)}]")
            output_blocks.append("\n".join(block_lines))

        # После последнего интервала
        if intervals and format_fraction(intervals[-1][1]) != "∞":
            output_blocks.append(f"({format_fraction(intervals[-1][1])}, ∞]\nрешений нет")

# в конце вместо formatted_output
        return {"intervals": [{"text": block} for block in output_blocks]}


    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))
