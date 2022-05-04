use concrete::*;
fn matrix_multi(matrix: &Vec<Vec<f64>>, vector: &Vec<f64>) -> Vec<f64> {
    let mut result: Vec<f64> = Vec::new();
    for i in 0..matrix.len() {
        let row: &Vec<f64> = &matrix[i];
        let mut v: f64 = 0.0;
        for j in 0..row.len() {
            let m: f64 = row[j];
            v += m * vector[j];
        }
        result.push(v);
    }
    return result;
}
fn invers_number_using_newton_raphson(a: f64) -> f64 {
    // the initial value should be silghtly greater than -1/((d+1)*(n + 1) * max * max)
    // where max is the max value in the X matrix
    // this is due how temp is calculated
    // a is negative so we should chose starting value to be slightly > 1/a
    // otherwise we will not converge into a value, because we are closer to critical point
    // https://blogs.sas.com/content/iml/2015/06/24/sensitivity-newtons-method.html#:~:text=If%20you%20provide%20a%20guess,away%20from%20the%20initial%20guess. 

    let mut xk: f64 = -0.001;
    for _ in 0..10 {
        xk = xk * (2.0 - a * xk);
        //println!("xk: {}", xk);
    }
    //println!("inv of : {} is {}", a, xk);
    return xk;
}

fn vect_mult(x: &Vec<f64>, y: &Vec<f64>) -> f64 {
    let mut r: f64 = 0.0;
    for i in 0..x.len() {
        r += x[i] * y[i];
    }
    return r;
}
fn vect_add(x: &Vec<f64>, y: &Vec<f64>) -> Vec<f64> {
    let mut r: Vec<f64> = Vec::new();
    for i in 0..x.len() {
        r.push(x[i] + y[i]);
    }
    return r;
}
fn vect_mult_scalar(x: &Vec<f64>, y: f64) -> Vec<f64> {
    let mut r: Vec<f64> = Vec::new();
    for i in 0..x.len() {
        r.push(x[i] * y);
    }
    return r;
}
fn main() -> Result<(), CryptoAPIError> {
    let n = 5;
    let d = 3;
    //A BAD APPROIMATION IN THE HESSIAN INVERSE MATRIX
    //WILL CAUSE US TO HAVE A SLOW PROGRESS TO THE OPTIMUM SOLUTION FOR BETA
    //THIS IS WHY WE CAN TOLERATE TO DO MANY ITERATIONS IN CALCULATING THE INVERSE MATRIX
    //WE WILL GAIN IN CALCULATION OF THE BETA
    let nbr_iters = 100;
    let x: Vec<Vec<f64>> = vec![
        vec![0.0, 0.0, 4.0, 4.0],
        vec![0.0, 0.0, 4.0, 4.0],
        vec![0.0, 0.0, 4.0, 4.0],
        vec![1.0, 1.0, 0.0, 0.0],
        vec![1.0, 1.0, 0.0, 0.0],
    ];

    let y: Vec<f64> = vec![1.0, 1.0, 1.0, -1.0, -1.0];
    let mut h_tild: Vec<Vec<f64>> = vec![vec![0.0; d + 1]; d + 1];
    let mut h_tild_inv: Vec<Vec<f64>> = vec![vec![0.0; d + 1]; d + 1];
    let mut beta: Vec<f64> = vec![0.001; d + 1];
    let mut sum: Vec<f64> = vec![0.0; n];
    for i in 0..n {
        for j in 0..d + 1 {
            sum[i] = sum[i] + x[i][j];
        }
    }

    for j in 0..d + 1 {
        let mut temp = 0.0;
        for i in 1..n {
            temp = temp + x[i][j] * sum[i];
        }
        h_tild[j][j] = -temp / 4.0;
        h_tild_inv[j][j] = invers_number_using_newton_raphson(h_tild[j][j]);
    }

    let mut deltas_history = Vec::new();
    for _ in 1..nbr_iters {
        let mut g: Vec<f64> = vec![0.0; d + 1];
        for i in 1..n {
            let a = (0.5 - 0.25 * y[i] * vect_mult(&x[i], &beta)) * y[i];
            g = vect_add(&g, &vect_mult_scalar(&x[i], a));
        }

        let delta = matrix_multi(&h_tild_inv, &g);
        beta = vect_add(&beta, &vect_mult_scalar(&delta, -1.0));
        deltas_history.push(delta);
    }

    println!("beta: {:?}", beta);
    Ok(())
}
