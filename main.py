from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from v4 import ParametricSimplex
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
    b_double: list[float]          # b'' (b_double)

def convert_to_fractions(data):
    """Конвертирует float в Fraction для точных вычислений."""
    return [
        [Fraction(str(x)) for x in row]
        if isinstance(row, list)
        else Fraction(str(row))
        for row in data
    ]

def format_fraction(value):
    if value == -math.inf:
        return "-∞"
    if value == math.inf:
        return "∞"
    if isinstance(value, Fraction):
        return float(value.numerator / value.denominator)
    return float(value)


@app.post("/solve-assignment")
async def solve_problem(data: ProblemData):
    try:
        print("Received data:", data.dict())
        # Извлечение и конвертация данных
        A = convert_to_fractions(data.matrix)
        c = convert_to_fractions(data.cost)
        b_prime = convert_to_fractions(data.b)
        b_double = convert_to_fractions(data.b_double)

        

        # Решение задачи
        ps = ParametricSimplex(A, b_prime, b_double, c)
        intervals = ps.solve()

        results = []
        for interval in intervals:
            L, U, alpha, beta, basis = interval
            
            # Форматирование lambda-диапазона (числовой формат)
            L_num = float(L) if L != -math.inf else "-inf"
            U_num = float(U) if U != math.inf else "inf"

            # Вычисление коэффициентов целевой функции
            c0 = sum(c[var] * alpha[i] for i, var in enumerate(basis) if var < len(c))
            c1 = sum(c[var] * beta[i] for i, var in enumerate(basis) if var < len(c))

            # Фильтрация ненулевых переменных
            X_expr = {}
            for j in range(len(c)):
                if j in basis:
                    idx = basis.index(j)
                    a_j = alpha[idx]
                    b_j = beta[idx]
                else:
                    a_j = Fraction(0)
                    b_j = Fraction(0)
                
                if a_j != 0 or b_j != 0:
                    a_str = format_fraction(a_j)
                    b_str = format_fraction(b_j)
                    expr = f"{a_str} + {b_str}λ" if b_j >= 0 else f"{a_str} - {format_fraction(abs(b_j))}λ"
                    X_expr[f"x{j+1}"] = expr

            results.append({
                "lambda_range": {"lower": L_num, "upper": U_num},
                "objective_coefficients": {
                    "c0": format_fraction(c0),
                    "c1": format_fraction(c1)
                },
                "non_zero_variables": X_expr
            })

        return {"status": "success", "result": results}

    except Exception as e:
        print("Ошибка на сервере:", str(e))  # Вывод полного трейса
        raise HTTPException(status_code=400, detail=str(e))