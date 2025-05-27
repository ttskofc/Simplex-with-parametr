document.addEventListener("DOMContentLoaded", function () {
    const initialRows = 3;
    const initialColumns = 5;
    createCost(initialColumns);
    createMatrix(initialRows, initialColumns);
    createB(initialRows);
    createBshtr(initialRows)
    setupKeyboardNavigation();

    document
        .getElementById("update-size")
        .addEventListener("click", updateMatrixSize);
    document
        .getElementById("submit-matrix")
        .addEventListener("click", submitMatrix);

    document.getElementById("matrix-size-rows").addEventListener("keydown", function (e) {
        if (e.key === "Enter") {
            updateMatrixSize();
        }
    });
    document.getElementById("matrix-size-columns").addEventListener("keydown", function (e) {
        if (e.key === "Enter") {
            updateMatrixSize();
        }
    });
});

// Создание матриц

function createMatrix(rows, columns) {
    const table = document.getElementById("matrix-table");
    table.innerHTML = "";

    for (let i = 0; i < rows; i++) {
        const row = document.createElement("tr");

        for (let j = 0; j < columns; j++) {
            const cell = document.createElement("td");
            const input = document.createElement("input");
            input.type = "number";
            input.className = "matrix-input";
            input.value = "0";
            input.dataset.row = i;
            input.dataset.col = j;
            cell.appendChild(input);
            row.appendChild(cell);
        }

        table.appendChild(row);
    }
}

function createCost(size) {
    const table = document.getElementById("cost-table");
    table.innerHTML = "";

    const row = document.createElement("tr");

    for (let i = 0; i < size; i++) {
        const cell = document.createElement("td");
        const input = document.createElement("input");
        input.type = "number";
        input.className = "cost-input";
        input.value = "0";
        input.dataset.col = i;
        cell.appendChild(input);
        row.appendChild(cell);
    }

    table.appendChild(row);

    const firstInput = document.querySelector(".cost-input");
    if (firstInput) {
        firstInput.focus();
    }
}

function createB(size) {
    const table = document.getElementById("b-table");
    table.innerHTML = "";

    for (let i = 0; i < size; i++) {
        const row = document.createElement("tr");

        const cell = document.createElement("td");
        const input = document.createElement("input");
        input.type = "number";
        input.className = "b-input";
        input.value = "0";
        input.dataset.row = i;
        cell.appendChild(input);
        row.appendChild(cell);

        table.appendChild(row);
    }
}

function createBshtr(size) {
    const table = document.getElementById("b-shtr-table");
    table.innerHTML = "";

    for (let i = 0; i < size; i++) {
        const row = document.createElement("tr");

        const cell = document.createElement("td");
        const input = document.createElement("input");
        input.type = "number";
        input.className = "b-shtr-input";
        input.value = "0";
        input.dataset.row = i;
        // input.dataset.col = j;
        cell.appendChild(input);
        row.appendChild(cell);

        table.appendChild(row);
    }
}

// Стрелки
function setupKeyboardNavigation() {
    const tableIds = ['matrix-table', 'cost-table', 'b-table', 'b-shtr-table'];
    tableIds.forEach(id => {
        const table = document.getElementById(id);
        if (table) {
            table.addEventListener('keydown', function(e) {
                if (['ArrowUp', 'ArrowDown', 'ArrowLeft', 'ArrowRight'].includes(e.key)) {
                    handleArrowNavigation(e);
                }
            });
        }
    });
}

// Обновленная функция обработки навигации
function handleArrowNavigation(e) {
    e.preventDefault();
    const currentInput = e.target;
    
    if (!currentInput.matches('.matrix-input, .cost-input, .b-input, .b-shtr-input')) return;
    
    const table = currentInput.closest('table');
    const tableId = table.id;
    
    let row, col;
    switch (tableId) {
        case 'cost-table':
            row = 0;
            col = parseInt(currentInput.dataset.col);
            break;
        case 'b-table':
        case 'b-shtr-table':
            row = parseInt(currentInput.dataset.row);
            col = 0;
            break;
        case 'matrix-table':
            row = parseInt(currentInput.dataset.row);
            col = parseInt(currentInput.dataset.col);
            break;
        default:
            return;
    }
    
    let rows, columns;
    switch (tableId) {
        case 'cost-table':
            rows = 1;
            columns = table.querySelectorAll('td').length;
            break;
        case 'matrix-table':
            rows = table.querySelectorAll('tr').length;
            columns = rows > 0 ? table.querySelector('tr').cells.length : 0;
            break;
        case 'b-table':
        case 'b-shtr-table':
            rows = table.querySelectorAll('tr').length;
            columns = 1;
            break;
        default:
            return;
    }
    
    let nextRow = row;
    let nextCol = col;
    
    switch (e.key) {
        case 'ArrowUp':
            if (tableId === 'cost-table') return;
            nextRow = Math.max(0, row - 1);
            break;
        case 'ArrowDown':
            if (tableId === 'cost-table') return;
            nextRow = Math.min(rows - 1, row + 1);
            break;
        case 'ArrowLeft':
            if (tableId === 'b-table' || tableId === 'b-shtr-table') return;
            nextCol = Math.max(0, col - 1);
            break;
        case 'ArrowRight':
            if (tableId === 'b-table' || tableId === 'b-shtr-table') return;
            nextCol = Math.min(columns - 1, col + 1);
            break;
        default:
            return;
    }
    
    let nextInput;
    switch (tableId) {
        case 'cost-table':
            nextInput = table.querySelector(`.cost-input[data-col="${nextCol}"]`);
            break;
        case 'b-table':
            nextInput = table.querySelector(`.b-input[data-row="${nextRow}"]`);
            break;
        case 'b-shtr-table':
            nextInput = table.querySelector(`.b-shtr-input[data-row="${nextRow}"]`);
            break;
        case 'matrix-table':
            nextInput = table.querySelector(`.matrix-input[data-row="${nextRow}"][data-col="${nextCol}"]`);
            break;
    }
    
    if (nextInput) {
        nextInput.focus();
        nextInput.select();
    }
}


// Обновлялка размера
function updateMatrixSize() {
    const newRow = parseInt(document.getElementById("matrix-size-rows").value);
    const newColumn = parseInt(document.getElementById("matrix-size-columns").value);
    
    createCost(newColumn);
    createMatrix(newRow, newColumn);
    createB(newRow)
    createBshtr(newRow)
}


async function submitMatrix() {
    const data = getMatrixData();
    
    
    try {
        const response = await fetch("http://localhost:8000/solve-assignment", {
            method: "POST",
            headers: {
                "Content-Type": "application/json",
            },
            body: JSON.stringify(data),
        });

        const result = await response.json();
        
        if (result.status === "success") {
            displayResult(result.result);
        } else {
            throw new Error(result.message || "Unknown error");
        }
    } catch (error) {
        console.error("Ошибка:", error);
        const outputDiv = document.getElementById("output");
        outputDiv.innerHTML = `<div class="error-message">Ошибка: ${error.message}</div>`;
    }
}



function displayResult(result) {
    const outputDiv = document.getElementById("output");
    let html = '';
    
    result.forEach(interval => {
        html += `
            <div class="interval">
                <h3>λ ∈ ${interval.lambda_range}</h3>
                <p>F(λ) = ${interval.objective}</p>
                <div class="variables">
                    <span>X = [</span>
                    ${interval.variables.map(v => `
                        <span class="variable">${v}</span>
                    `).join(', ')}
                    <span>]</span>
                </div>
            </div>
        `;
    });
    
    outputDiv.innerHTML = html;
}


// данные матрицы
function getMatrixData() {
    const row = document.querySelectorAll("#matrix-table tr").length;
    const column = document.querySelectorAll("#matrix-table tr:first-child td").length;

    const matrix = [];
    const cost = [];
    const b = [];
    const bshtr = []; 

    for (let i = 0; i < row; i++) {
        matrix[i] = [];
        for (let j = 0; j < column; j++) {
            const input = document.querySelector(
                `.matrix-input[data-row="${i}"][data-col="${j}"]`
            );
            matrix[i][j] = parseFloat(input.value) || 0;
        }
    }

    for (let j = 0; j < column; j++) {
        const input = document.querySelector(
            `.cost-input[data-col="${j}"]`
        );
        cost[j] = parseFloat(input.value) || 0;
    }

    for (let i = 0; i < row; i++) {
        const input = document.querySelector(
            `.b-input[data-row="${i}"]` 
        );
        b[i] = parseFloat(input.value) || 0; 
    }

    for (let i = 0; i < row; i++) {
        const input = document.querySelector(
            `.b-shtr-input[data-row="${i}"]` 
        );
        bshtr[i] = parseFloat(input.value) || 0; 
    }

    const data = {
        matrix: matrix,
        cost: cost,
        b: b,
        b_double: bshtr, 
    };

    // alert(JSON.stringify(data, null, 2));
    return data;
}