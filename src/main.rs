use concrete::*;
//use ndarray::arr2;
//use ndarray::Array;
//use ndarray::Array1;
/**
 * 
n=10 ; %number of training records
d = 10; %number of features
nbr_iters = 10;
x = zeros(n, d + 1);
y = zeros(n);


h_tild = zeros(d+1, d+1);
h_tild_inv = zeros(d+1, d+1);

beta = 0.001 * ones(d + 1, 1);
sum = zeros(n, 1);
for i = 1 : n 
    for j = 1 : d + 1
        sum(i) = sum(i) + x(i, j);
    end
end
for j=1 : d+1
    temp=0;
    for i = 1 : n
    temp = temp + x(i, j) * sum(i);
    end
    h_tild(j,j) = 1/4 * temp;
    h_tild_inv(j,j) = 2*u-h_tild(j,j) *power(u,2);
end
deltasHistory = zeros(nbr_iters, 1);
for k = 1 : nbr_iters
    for i = 1 : n
		g= g + (1/2-1/4 * y (i) * x(i,:)* beta)*y(i)*x(i,:);
    end 
    delta = h_tild_inv * g;
    deltasHistory(k) = delta;
	beta = beta - delta;
end
plot(deltasHistory);
*/
fn matrix_multi(matrix: &Vec<Vec<f32>>, vector: &Vec<f32>) -> Vec<f32>{
    let mut result :Vec<f32> = Vec::new();
    for i in 0.. matrix.len(){
        let row: &Vec<f32> = &matrix[i];
        let mut v : f32 = 0.0;
        for j in 0.. row.len(){
            let m : f32= row[j];
            v += m * vector[j];
        }
        result.push(v);
    }
    return result;
}
fn invers_number_using_newton_raphson(a:f32)->f32{
    //let mut xk:f32 = -0.00001;
    //for _ in 0..20{
    //   xk = xk * (2.0 - a*xk);
    //   println!("xk: {}", xk);
    //}
    //println!("inv of : {} is {}", a, xk);
    //return xk;
    return 1.0/a;
}

fn vect_mult(x: &Vec<f32>, y: &Vec<f32>)->f32{
    let mut r:f32 = 0.0;
    for i in 0..x.len(){
        r += x[i] * y[i];
    }
    return r;
}
fn vect_add(x: &Vec<f32>, y: &Vec<f32>)->Vec<f32>{
    let mut r:Vec<f32> = Vec::new();
    for i in 0..x.len(){
        r.push(x[i] + y[i]);
    }
    return r;
}
fn vect_mult_scalar(x: &Vec<f32>, y: f32)->Vec<f32>{
    let mut r:Vec<f32> = Vec::new();
    for i in 0..x.len(){
        r.push(x[i] * y);
    }
    return r;
}
fn main() -> Result<(), CryptoAPIError> {

    let n = 5;
    let d = 3;
    let nbr_iters = 100;
    let x: Vec<Vec<f32>>= vec![
                                vec![0.0,0.0,4.0,4.0],
                                vec![0.0,0.0,4.0,4.0],
                                vec![0.0,0.0,4.0,4.0],
                                vec![1.0,1.0,0.0,0.0],
                                vec![1.0,1.0,0.0,0.0],
    ];
    let y: Vec<f32>= vec![1.0,1.0,1.0,-1.0,-1.0];
    let mut h_tild: Vec<Vec<f32>>= vec![vec![0.0; d+1]; d+1];
    let mut h_tild_inv: Vec<Vec<f32>>= vec![vec![0.0; d+1]; d+1];
    let mut beta: Vec<f32>= vec![0.001; d+1];
    let mut sum: Vec<f32>= vec![0.0; n];
    for i in 0..n{
        for j in 0..d+1{
            sum[i] = sum[i] + x[i][j];
        }
    }
    for j in 0..d+1{
        let mut temp=0.0;
        for i in 1..n{
            temp = temp + x[i][j] * sum[i];
        }
        h_tild[j][j] =  - temp / 4.0;
        h_tild_inv[j][j] = invers_number_using_newton_raphson(h_tild[j][j]);
    }

    let mut deltas_history = Vec::new();
    for _ in 1 .. nbr_iters{
        let mut g : Vec<f32>= vec![0.0; d+1];
        //println!("==================================");
        for i in 1 ..n {
            //println!("multi x, beta: {:?}", vect_mult(&x[i], &beta));
            let a = (0.5-0.25 * y[i] * vect_mult(&x[i], &beta))*y[i];
            //println!("g: {:?}", g);
            g = vect_add(&g,  &vect_mult_scalar(&x[i], a));
        }
        //MAtrix multi
        let delta = matrix_multi(&h_tild_inv , &g);
        beta = vect_add(&beta , &vect_mult_scalar (&delta, -1.0));
        deltas_history.push(delta);
    }

    println!("beta: {:?}", beta);            // prints [6.26]
    //println!("h_tild_inv: {:?}", h_tild_inv);            // prints [6.26]
    //println!("delta history : {:?}", deltas_history);            // prints [6.26]

    Ok(())

}
