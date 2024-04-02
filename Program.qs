namespace QIEAnet {
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Convert;
    @EntryPoint()
    operation Main() : Double {
        let n = 10; //кол-во вершин
        let k = 2; // 2 в степени k >= корень из n (2 в k - (max) кол-во сообществ)
        let r = 12; // Кол-во ребер
        // Матрица смежности графа
        let a = [[0, 0, 0, 0, 1, 0, 0, 0, 1, 1],
                 [0, 0, 1, 1, 0, 0, 0, 0, 0, 1],
                 [0, 1, 0, 0, 0, 0, 0, 0, 0, 1],
                 [0, 1, 0, 0, 0, 0, 0, 0, 0, 1],
                 [1, 0, 0, 0, 0, 1, 0, 0, 1, 0],
                 [0, 0, 0, 0, 1, 0, 1, 1, 0, 0],
                 [0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                 [1, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                 [1, 1, 1, 1, 0, 0, 0, 0, 0, 0]];
        mutable t = 0; // Текущая эволюция
        let q = 60; // Угол поворота при эволюции
        let m = n * n; // Кол-во кв индивидов
        let maxGen = n * n; // Кол-во эволюций
        // Создание начальных квантовых состояний для популяции
        mutable Q = CreateQ(n, m, k);
        mutable P = CreateP(n, m, k, Q, q);
        mutable C = CreateC(n, m, k, P);
        mutable F = CreateF(a, r, n, m, k, C);
        mutable Best = FindBest(F);
        mutable BestC = FindBestC(F, C, Best);
        mutable BestP = FindBestP(F, C, P, Best);
        // Основной цикл эволюции
        while (t <= maxGen){
            set t += 1;
            // Генерация новых решений
            mutable P = CreateP(n, m, k, Q, q);
            mutable C = CreateC(n, m, k, P);
            // Оценка новых решений
            mutable F = CreateF(a, r, n, m, k, C);
            mutable newBest = FindBest(F);
            mutable newBestC = FindBestC(F, C, Best);
            mutable newBestP = FindBestP(F, C, P, Best);
            // Обновление решения если найдено лучше
            if (newBest > Best) {
                set Best = newBest;
                set BestC = newBestC;
                set BestP = newBestP;
                // Эволюция квантовых индивидов
                set Q = UpdateUsingGates(Q, P, BestP, n, m, k);
            }

        }
        // Вывод лучшего разделения на сообщества
        Message($"{BestC}");
        // Возвращение значения лучшей модульности
        return Best;
    }
    operation CreateP(n : Int, m : Int, k : Int, Q : Int[][], q : Int) : Int[][] {
        mutable P = [];
        for i in 0 .. m - 1 {
            use newRow = Qubit[n * k];
            mutable P_row = [];
            for j in 0 .. n * k - 1 {
                H(newRow[j]);
                for j in 0 .. n * k - 1 {
                    mutable rotate = 3.142 / (180.0 / IntAsDouble(q)) * IntAsDouble(Q[i][j]);
                    if (rotate > 3.142) {
                        set rotate = 3.142;
                    }
                    if (rotate < -3.142) {
                        set rotate = -3.142;
                    }
                    Ry(rotate, newRow[j]);
                }
                let res = M(newRow[j]);
                if (res == One){
                    set P_row += [1];
                } else {
                    set P_row += [0];
                }
                Reset(newRow[j]);
            }
            ResetAll(newRow);
            set P += [P_row];
        }
        return P;
    }
    operation CreateC(n : Int, m : Int, k : Int, P : Int[][]) : Int[][] {
        mutable C = [];
        for i in 0 .. m - 1 {
            mutable C_row = [];
            for j in 0 .. n - 1{
                set C_row += [IntFromDouble(P[i][j * k .. (j * k + k - 1)])];
            }
            set C += [C_row];
        }
        return C;
    }
    operation CreateF(a : Int[][], r : Int, n : Int, m : Int, k : Int, C : Int[][]) : Double[] {
        mutable F = [];
        mutable F_row = 0.0;
         for i in 0 .. m - 1 {
            for j in 0 .. n - 1{
                set F_row = ModuleFunk(C[i], a, r);
            }
            set F += [F_row];
        }
        return F;
    }
    operation CreateQ(n : Int, m : Int, k : Int) : Int[][] {
        mutable Q = [];
         for i in 0 .. m - 1 {
            mutable Q_row = [];
            for j in 0 .. n * k - 1{
                set Q_row += [0];
            }
            set Q += [Q_row];
        }
        return Q;
    }
    operation ModuleFunk(array1 : Int[], array2 : Int[][], m : Int) : Double {
        mutable answ = 0.0;
        mutable a = 0.0;
        mutable k = [];
        for i in 0 .. Length(array2) - 1 {
            mutable t = 0.0;
            for j in 0 .. Length(array2) - 1 {
                if (array2[i][j] == 1)
                {
                    set t += 1.0;
                }
            }
            set k += [t];
        }
        for i in 0 .. Length(array2) - 1{
            for j in 0 .. Length(array2) - 1 {
                if (i != j){
                    if (array2[i][j] == 1) {
                        set a = 1.0;
                    } else {
                        set a = 0.0;
                    }
                    if (array1[i] == array1[j]){
                        set answ += a - ((k[i] * k[j]) / (2.0 * IntAsDouble(m)));
                    }
                }
            }
        }
        return answ / (2.0 * IntAsDouble(m));
    }
    operation UpdateUsingGates(Q : Int[][], P : Int[][], BestP : Int[], n : Int, m : Int, k : Int) : Int[][] {
        mutable newQ = [];
        for i in 0 .. m - 1 {
            mutable newQ_row = [];
            for j in 0 .. n - 1{
                for t in 0 .. k - 1{
                    if (BestP[j * k + t] != P[i][j * k + t]) {
                        set newQ_row += [0];
                    } else {
                        if (P[i][j * k + t] == 0) {
                            set newQ_row += [1];
                        } else {
                            set newQ_row += [-1];
                        }
                    }
                }
            }
            set newQ += [newQ_row];
        }
        return newQ;
    }
    operation FindBest(array : Double[]) : Double {
        mutable answ = -10000.0;
        for i in 0 .. Length(array) - 1 {
            if (answ < array[i]) {
                set answ = array[i];
            }
        }
        return answ;
    }
    operation FindBestC(array : Double[], array2 : Int[][], max : Double) : Int[] {
        mutable ar = array2[0];
        for i in 0 .. Length(array) - 1 {
            if (max == array[i]) {
                return array2[i];
            }
        }
        return ar;
    }
    operation FindBestP(array : Double[], array2 : Int[][], array3 : Int[][], max : Double) : Int[] {
        mutable ar = array2[0];
        for i in 0 .. Length(array) - 1 {
            if (max == array[i]) {
                return array3[i];
            }
        }
        return array3[0];
    }
    operation IntFromDouble(array : Int[]) : Int {
        mutable answ = 0;
        mutable t = 0;
        for i in Length(array) - 1 .. -1 .. 0 {
            mutable k = 1;
            for j in 1 .. t {
                set k *= 2;
            }
            set answ += array[i] * k;
            set t += 1;
        }
        return answ;
    }
}
