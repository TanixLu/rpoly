const FLT_MIN: f64 = f32::MIN_POSITIVE as f64;
const FLT_MAX: f64 = f32::MAX as f64;
const DBL_EPSILON: f64 = f64::EPSILON;

fn fabs(x: f64) -> f64 {
    x.abs()
}

fn sqrt(x: f64) -> f64 {
    x.sqrt()
}

fn log(x: f64) -> f64 {
    x.ln()
}

fn exp(x: f64) -> f64 {
    x.exp()
}

fn pow(x: f64, p: i32) -> f64 {
    x.powi(p)
}

fn cos(x: f64) -> f64 {
    x.cos()
}

fn sin(x: f64) -> f64 {
    x.sin()
}

fn rpoly_ak1<const MDP1: usize>(
    op: &[f64; MDP1],
    Degree: &mut usize,
    zeror: &mut [f64; MDP1],
    zeroi: &mut [f64; MDP1],
) {
    const RADFAC: f64 = std::f64::consts::PI / 180.0;
    const lb2: f64 = std::f64::consts::LN_2;
    const lo: f64 = FLT_MIN / DBL_EPSILON;
    let cosr = cos(94.0 * RADFAC);
    let sinr = sin(94.0 * RADFAC);

    // Do a quick check to see if leading coefficient is 0
    if op[0] != 0.0 {
        let mut N = *Degree;
        let mut xx = std::f64::consts::FRAC_1_SQRT_2; // = 0.70710678
        let mut yy = -xx;

        // Remove zeros at the origin, if any
        let mut j = 0;
        while op[N] == 0.0 {
            zeroi[j] = 0.0;
            zeror[j] = 0.0;
            N -= 1;
            j += 1
        }

        let mut NN = N + 1;

        // Make a copy of the coefficients
        let mut p = [0.0; MDP1];
        for i in 0..NN {
            p[i] = op[i]
        }

        while N >= 1 {
            // Main loop
            // Start the algorithm for one zero
            if N <= 2 {
                // Calculate the final zero or pair of zeros
                if N < 2 {
                    zeror[*Degree - 1] = -(p[1] / p[0]);
                    zeroi[*Degree - 1] = 0.0;
                } else {
                    let [mut sr, mut si, mut lr, mut li] = [0.0; 4];
                    Quad_ak1(p[0], p[1], p[2], &mut sr, &mut si, &mut lr, &mut li);
                    zeror[*Degree - 2] = sr;
                    zeroi[*Degree - 2] = si;
                    zeror[*Degree - 1] = lr;
                    zeroi[*Degree - 1] = li;
                }
                break;
            }

            // Find the largest and smallest moduli of the coefficients

            let mut moduli_max = 0.0;
            let mut moduli_min = FLT_MAX;

            for i in 0..NN {
                let x = fabs(p[i]);
                if x > moduli_max {
                    moduli_max = x;
                }
                if x != 0.0 && x < moduli_min {
                    moduli_min = x;
                }
            }

            // Scale if there are large or very small coefficients
            // Computes a scale factor to multiply the coefficients of the polynomial. The scaling
            // is done to avoid overflow and to avoid undetected underflow interfering with the
            // convergence criterion.
            // The factor is a power of the base.

            let mut sc = lo / moduli_min;

            if ((sc <= 1.0) && (moduli_max >= 10.0)) || ((sc > 1.0) && (FLT_MAX / sc >= moduli_max))
            {
                if sc == 0.0 {
                    sc = FLT_MIN;
                    let l = (log(sc) / lb2 + 0.5) as i32;
                    let factor = pow(2.0, l);
                    if factor != 1.0 {
                        for i in 0..NN {
                            p[i] *= factor;
                        }
                    }
                }
            }

            // Compute lower bound on moduli of zeros

            let mut pt = [0.0; MDP1];
            for i in 0..NN {
                pt[i] = fabs(p[i]);
            }
            pt[N] = -(pt[N]);

            let NM1 = N - 1;

            // Compute upper estimate of bound
            let mut x = exp((log(-pt[N]) - log(pt[0])) / N as f64);

            if pt[NM1] != 0.0 {
                // If Newton step at the origin is better, use it
                let xm = -pt[N] / pt[NM1];
                if xm < x {
                    x = xm;
                }
            }

            // Chop the interval (0, x) until ff <= 0

            let mut xm = x;
            loop {
                x = xm;
                xm = 0.1 * x;
                let mut ff = pt[0];
                for i in 1..NN {
                    ff = ff * xm + pt[i];
                }
                if ff <= 0.0 {
                    break;
                }
            }

            let mut dx = x;

            // Do Newton iteration until x converges to two decimal places
            while fabs(dx / x) > 0.005 {
                let mut ff = pt[0];
                let mut df = pt[0];
                for i in 1..N {
                    ff = x * ff + pt[i];
                    df = x * df + ff;
                }
                ff = x * ff + pt[N];
                dx = ff / df;
                x -= dx;
            }

            let bnd = x;

            // Compute the derivative as the initial K polynomial and do 5 steps with no shift

            let mut K = [0.0; MDP1];
            for i in 1..N {
                K[i] = ((N - i) as f64) * p[i] / (N as f64);
            }
            K[0] = p[0];

            let mut aa = p[N];
            let mut bb = p[NM1];
            let mut zerok = K[NM1] == 0.0;

            for jj in 0..5 {
                let cc = K[NM1];
                if zerok {
                    // Use unscaled form of recurrence
                    for i in 0..NM1 {
                        j = NM1 - i;
                        K[j] = K[j - 1];
                    }
                    K[0] = 0.0;
                    zerok = K[NM1] == 0.0;
                } else {
                    // Used scaled form of recurrence if value of K at 0 is nonzero
                    let t = -aa / cc;
                    for i in 0..NM1 {
                        j = NM1 - i;
                        K[j] = t * K[j - 1] + p[j];
                    }
                    K[0] = p[0];
                    zerok = fabs(K[NM1]) <= fabs(bb) * DBL_EPSILON * 10.0;
                }
            }

            // Save K for restarts with new shifts
            let mut temp = [0.0; MDP1];
            for i in 0..N {
                temp[i] = K[i];
            }

            // Loop to select the quadratic corresponding to each new shift

            let mut jj = 0;

            for jj in 1..=20 {
                // Quadratic corresponds to a double shift to a non-real point and its
                // complex conjugate. The point has modulus BND and amplitude rotated
                // by 94 degrees from the previous shift.

                let xxx = -(sinr * yy) + cosr * xx;
                yy = sinr * xx * cosr * yy;
                xx = xxx;
                let sr = bnd * xx;

                // Second stage calculation, fixed quadratic

                let mut NZ = 0;
                let mut qp = [0.0; MDP1];
                let [mut lzi, mut lzr, mut szi, mut szr] = [0.0; 4];
                Fxshfr_ak1(
                    20 * jj,
                    &mut NZ,
                    sr,
                    bnd,
                    &mut K,
                    N,
                    &mut p,
                    NN,
                    &mut qp,
                    &mut lzi,
                    &mut lzr,
                    &mut szi,
                    &mut szr,
                );

                if NZ != 0 {
                    // The second stage jumps directly to one of the third stage iterations and
                    // returns here if successful. Deflate the polynomial, store the zero or
                    // zeros, and return to the main algorithm.

                    j = *Degree - N;
                    zeror[j] = szr;
                    zeroi[j] = szi;
                    NN = NN - NZ as usize;
                    N = NN - 1;

                    for i in 0..NN {
                        p[i] = qp[i];
                    }
                    if NZ != 1 {
                        zeror[j + 1] = lzr;
                        zeroi[j + 1] = lzi;
                    }
                    break;
                } else {
                    // If the iteration is unsuccessful, another quadratic is chosen after restoring K
                    for i in 0..N {
                        K[i] = temp[i];
                    }
                }
            } // End for jj

            // Return with failure if no convergence with 20 shifts
            if jj > 20 {
                // TODO: replace panic with Error
                panic!("Failure. No convergence after 20 shifts. Program terminated.");
            }
        } // End while (N >= 1)
    } else {
        // *Degree = 0;
        // TODO: remove this panic
        panic!("The leading coefficient is zero. No further action taken. Program terminated.");
    }
}

fn Fxshfr_ak1<const MDP1: usize>(
    L2: i32,
    NZ: &mut i32,
    sr: f64,
    bnd: f64,
    K: &mut [f64; MDP1],
    N: usize,
    p: &mut [f64; MDP1],
    NN: usize,
    qp: &mut [f64; MDP1],
    lzi: &mut f64,
    lzr: &mut f64,
    szi: &mut f64,
    szr: &mut f64,
) {
    // Computes up to L2 fixed shift K-polynomials, testing for convergence in the linear or
    // quadratic case. Initiates one of the variable shift iterations and returns with the
    // number of zeros found.

    // L2 limit of fixed shift steps
    // NZ number of zeros found

    *NZ = 0;
    let mut betas = 0.25;
    let mut betav = 0.25;
    let mut u = -(2.0 * sr);
    let mut oss = sr;
    let mut v = bnd;
    let mut ovv = bnd;

    // Evaluate polynomial by synthetic division
    let mut a = 0.0;
    let mut b = 0.0;
    QuadSD_ak1(NN, u, v, p, qp, &mut a, &mut b);

    let [mut a1, mut a3, mut a7, mut c, mut d, mut e, mut f, mut g, mut h] = [0.0; 9];
    let mut qk = [0.0; MDP1];
    let mut tFlag = calcSC_ak1(
        N, a, b, &mut a1, &mut a3, &mut a7, &mut c, &mut d, &mut e, &mut f, &mut g, &mut h, K, u,
        v, &mut qk,
    );

    let [mut ovv, mut oss, mut otv, mut ots] = [0.0; 4];
    for j in 0..L2 {
        // Calculate next K polynomial and estimate v
        nextK_ak1(N, tFlag, a, b, a1, &mut a3, &mut a7, K, &mut qk, qp);
        tFlag = calcSC_ak1(
            N, a, b, &mut a1, &mut a3, &mut a7, &mut c, &mut d, &mut e, &mut f, &mut g, &mut h, K,
            u, v, &mut qk,
        );
        let mut ui = 0.0;
        let mut vi = 0.0;
        newest_ak1(
            tFlag, &mut ui, &mut vi, a, a1, a3, a7, b, c, d, f, g, h, u, v, K, N, p,
        );

        let mut vv = vi;

        // Estimate s

        let mut ss = if K[N - 1] != 0.0 {
            -(p[N] / K[N - 1])
        } else {
            0.0
        };

        let mut tv = 1.0;
        let mut ts = 1.0;

        if j != 0 && tFlag != 3 {
            // Compute relative measures of convergence of s and v sequences

            if vv != 0.0 {
                tv = fabs((vv - ovv) / vv);
            }
            if ss != 0.0 {
                ts = fabs((ss - oss) / ss);
            }

            // If decreasing, multiply the two most recent convergence measures

            let mut tvv = if tv < otv { tv * otv } else { 1.0 };
            let mut tss = if ts < ots { ts * ots } else { 1.0 };

            // Compare with convergence criteria

            let vpass = tvv < betav;
            let spass = tss < betas;

            if spass || vpass {
                // At least one sequence has passed the convergence test.
                // Store variables before iterating

                let mut svk = [0.0; MDP1];
                for i in 0..N {
                    svk[i] = K[i];
                }

                let mut s = ss;

                // Choose iteration according to the fastest converging sequence

                let mut vtry = false;
                let mut stry = false;
                let mut fflag = true;

                loop {
                    let mut iFlag = true; // Begin each loop by assuming RealIT will be called UNLESS iFlag changed below

                    if fflag && spass && (!vpass || tss < tvv) {
                        fflag = false;
                    } else {
                        if fflag {
                            fflag = false;
                        }

                        QuadIT_ak1(
                            N, NZ, ui, vi, szr, szi, lzr, lzi, qp, NN, &mut a, &mut b, p, &mut qk,
                            &mut a1, &mut a3, &mut a7, &mut d, &mut e, &mut f, &mut g, &mut h, K,
                        );

                        if *NZ > 0 {
                            return;
                        }

                        // Quadratic iteration has failed. Flag that it has been tried and decrease the
                        // convergence criterion

                        vtry = true;
                        betav *= 0.25;

                        // Try linear iteration if it has not been tried and the s sequence is converging
                        if stry || !spass {
                            iFlag = false;
                        } else {
                            for i in 0..N {
                                K[i] = svk[i];
                            }
                        }
                    }

                    if iFlag {
                        RealIT_ak1(&mut iFlag, NZ, &mut s, N, p, NN, qp, szr, szi, K, &mut qk);

                        if *NZ > 0 {
                            return;
                        }

                        // Linear iteration has failed. Flag that it has been tried and decrease the
                        // convergence criterion

                        stry = true;
                        betas *= 0.25;

                        if iFlag {
                            // If linear iteration signals an almost double real zero, attempt quadratic iteration

                            ui = -(s + s);
                            vi = s * s;
                            continue;
                        }
                    }

                    // Restore variables
                    for i in 0..N {
                        K[i] = svk[i];
                    }

                    // Try quadratic iteration if it has not been tried and the v sequence is converging

                    if !vpass || vtry {
                        break;
                    }
                } // End do-while loop

                // Re-compute qp and scalar values to continue the second stage
                QuadSD_ak1(NN, u, v, p, qp, &mut a, &mut b);
                tFlag = calcSC_ak1(
                    N, a, b, &mut a1, &mut a3, &mut a7, &mut c, &mut d, &mut e, &mut f, &mut g,
                    &mut h, K, u, v, &mut qk,
                );
            } // End if ((spass) || (vpass))
        } // End if ((j != 0) && (tFlag != 3))

        ovv = vv;
        oss = ss;
        otv = tv;
        ots = ts;
    } // End for j
}

fn RealIT_ak1(
    iFlag: &mut bool,
    NZ: &mut i32,
    sss: &mut f64,
    N: usize,
    p: &mut [f64],
    NN: usize,
    qp: &mut [f64],
    szr: &mut f64,
    szi: &mut f64,
    K: &mut [f64],
    qk: &mut [f64],
) {
    // Variable-shift H-polynomial iteration for a real zero

    // sss - starting iterate
    // NZ - number of zeros found
    // iFlag - flag to indicate a pair of zeros near real axis

    let mut i = 0;
    let mut j = 0;
    let nm1 = N - 1;

    let [mut ee, mut kv, mut mp, mut ms, mut omp, mut pv, mut s, mut t] = [0.0; 8];

    *NZ = 0;
    *iFlag = false;
    let mut s = *sss;

    loop {
        let mut pv = p[0];
        qp[0] = pv;

        // Evaluate p at s
        for i in 1..NN {
            pv = pv * s + p[i];
            qp[i] = pv;
        }

        let mut mp = fabs(pv);

        // Compute a rigorous bound on the error in evaluating p

        let mut ms = fabs(s);
        let mut ee = 0.5 * fabs(qp[0]);
        for i in 1..NN {
            ee = ee * ms + fabs(qp[i]);
        }

        // Iteration has converged sufficiently if the polynomial value is less than
        // 20 times this bound

        if mp <= 20.0 * DBL_EPSILON * (2.0 * ee - mp) {
            *NZ = 1;
            *szr = s;
            *szi = 0.0;
            break;
        }

        j += 1;

        // Stop iteration after 10 steps
        if j > 10 {
            break;
        }

        if j >= 2 {
            if (fabs(t) <= 0.001 * fabs(-t + s)) && (mp > omp) {
                // A cluster of zeros near the real axis has been encountered;
                // Return with iFlag set to initiate a quadratic iteration

                *iFlag = true;
                *sss = s;
                break;
            }
        }

        // Return if the polynomial value has increased significantly

        omp = mp;

        // Compute t, the next polynomial and the new iterate
        kv = K[0];
        qk[0] = kv;
        for i in 1..N {
            kv = kv * s + K[i];
            qk[i] = kv;
        }

        if fabs(kv) > fabs(K[nm1]) * 10.0 * DBL_EPSILON {
            // Use the scaled form of the recurrence if the value of K at s is non-zero
            t = -(pv / kv);
            K[0] = qp[0];
            for i in 1..N {
                K[i] = t * qk[i - 1] + qp[i];
            }
        } else {
            // Use unscaled form
            K[0] = 0.0;
            for i in 1..N {
                K[i] = qk[i - 2];
            }
        }

        kv = K[0];
        for i in 1..N {
            kv = kv * s + K[i];
        }

        t = if fabs(kv) > (fabs(K[nm1]) * 10.0 * DBL_EPSILON) {
            -(pv / kv)
        } else {
            0.0
        };

        s += t;
    }
}

fn QuadIT_ak1(
    N: usize,
    NZ: &mut i32,
    uu: f64,
    vv: f64,
    szr: &mut f64,
    szi: &mut f64,
    lzr: &mut f64,
    lzi: &mut f64,
    qp: &mut [f64],
    NN: usize,
    a: &mut f64,
    b: &mut f64,
    p: &mut [f64],
    qk: &mut [f64],
    a1: &mut f64,
    a3: &mut f64,
    a7: &mut f64,
    d: &mut f64,
    e: &mut f64,
    f: &mut f64,
    g: &mut f64,
    h: &mut f64,
    K: &mut [f64],
) {
}

fn newest_ak1(
    tFlag: i32,
    uu: &mut f64,
    vv: &mut f64,
    a: f64,
    a1: f64,
    a3: f64,
    a7: f64,
    b: f64,
    c: f64,
    d: f64,
    f: f64,
    g: f64,
    h: f64,
    u: f64,
    v: f64,
    K: &mut [f64],
    N: usize,
    p: &mut [f64],
) {
    *uu = 0.0;
    *vv = 0.0;

    if tFlag != 3 {
        let (a4, a5) = if tFlag != 2 {
            (a + u * b + h * f, c + (u + v * f) * d)
        } else {
            ((a + g) * f + h, (f + u) * c + v * d)
        };

        // Evaluate new quadratic coefficients
        let b1 = -K[N - 1] / p[N];
        let b2 = -(K[N - 2] + b1 * p[N - 1]) / p[N];
        let c1 = v * b2 * a1;
        let c2 = b1 * a7;
        let c3 = b1 * b1 * a3;
        let c4 = -(c2 + c3) + c1;
        let temp = -c4 + a5 + b1 * a4;
        if temp != 0.0 {
            *uu = -((u * (c3 + c2) + v * (b1 * a1 + b2 * a7)) / temp) + u;
            *vv = v * (1.0 + c4 / temp);
        }
    }
}

fn nextK_ak1(
    N: usize,
    tFlag: i32,
    a: f64,
    b: f64,
    a1: f64,
    a3: &mut f64,
    a7: &mut f64,
    K: &mut [f64],
    qk: &mut [f64],
    qp: &mut [f64],
) {
    // Computes the next K polynomials using the scalars computed in calcSC_ak1

    if tFlag == 3 {
        // Use unscaled form of the recurrence
        K[0] = 0.0;
        K[1] = 0.0;

        for i in 2..N {
            K[i] = qk[i - 2];
        }

        return;
    }

    let temp = if tFlag == 1 { b } else { a };

    if fabs(a1) > (10.0 * DBL_EPSILON * fabs(temp)) {
        // Use scaled form of the recurrence

        (*a7) /= a1;
        (*a3) /= a1;
        K[0] = qp[0];
        K[1] = -((*a7) * qp[0]) + qp[1];

        for i in 2..N {
            K[i] = -((*a7) * qp[i - 1]) + (*a3) * qk[i - 2] + qp[i];
        }
    } else {
        // If a1 is nearly zero, then use a special form of the recurrence

        K[0] = 0.0;
        K[1] = -(*a7) * qp[0];

        for i in 2..N {
            K[i] = -((*a7) * qp[i - 1]) + (*a3) * qk[i - 2];
        }
    }
}

fn calcSC_ak1(
    N: usize,
    a: f64,
    b: f64,
    a1: &mut f64,
    a3: &mut f64,
    a7: &mut f64,
    c: &mut f64,
    d: &mut f64,
    e: &mut f64,
    f: &mut f64,
    g: &mut f64,
    h: &mut f64,
    K: &mut [f64],
    u: f64,
    v: f64,
    qk: &mut [f64],
) -> i32 {
    // This routine calculates scalar quantities used to compute the next K polynomial and
    // new estimates of the quadratic coefficients.

    // calcSC - integer variable set here indicating how the calculations are normalized
    // to avoid overflow.

    let mut dumFlag = 3; // TYPE = 3 indicates the quadratic is almost a factor of K

    // Synthetic division of K by the quadratic 1, u, v
    QuadSD_ak1(N, u, v, K, qk, c, d);

    if fabs(*c) <= (100.0 * DBL_EPSILON * fabs(K[N - 1])) {
        if fabs(*d) <= (100.0 * DBL_EPSILON * fabs(K[N - 2])) {
            return dumFlag;
        }
    }

    *h = v * b;
    if fabs(*d) >= fabs(*c) {
        dumFlag = 2; // TYPE = 2 indicates that all formulas are divided by d
        *e = a / (*d);
        *f = (*c) / (*d);
        *g = u * b;
        *a3 = (*e) * ((*g) + a) + (*h) * (b / (*d));
        *a1 = -a + (*f) * b;
        *a7 = (*h) + ((*f) + u) * a;
    } else {
        dumFlag = 1; // TYPE = 1 indicates that all formulas are divided by c;
        *e = a / (*c);
        *f = (*d) / (*c);
        *g = (*e) * u;
        *a3 = (*e) * a + ((*g) + (*h) / (*c)) * b;
        *a1 = -(a * ((*d) / (*c))) + b;
        *a7 = (*g) * (*d) + (*h) * (*f) + a;
    }

    dumFlag
}

fn QuadSD_ak1(NN: usize, u: f64, v: f64, p: &mut [f64], q: &mut [f64], a: &mut f64, b: &mut f64) {
    // Divides p by the quadratic 1, u, v placing the quotient in q and the remainder in a, b

    *b = p[0];
    q[0] = *b;
    *a = -((*b) * u) + p[1];
    q[1] = *a;

    for i in 2..NN {
        q[i] = -((*a) * u + (*b) * v) + p[i];
        *b = *a;
        *a = q[i];
    }
}

fn Quad_ak1(a: f64, b1: f64, c: f64, sr: &mut f64, si: &mut f64, lr: &mut f64, li: &mut f64) {
    let mut b: f64 = 0.;
    let mut d: f64 = 0.;
    let mut e: f64 = 0.;
    *li = 0.0f64;
    *lr = *li;
    *si = *lr;
    *sr = *si;
    if a == 0.0 {
        *sr = if b1 != 0.0 { -(c / b1) } else { *sr };
        return;
    }
    if c == 0.0 {
        *lr = -(b1 / a);
        return;
    }
    b = b1 / 2.0f64;
    if fabs(b) < fabs(c) {
        e = if c >= 0.0 { a } else { -a };
        e = -e + b * (b / fabs(c));
        d = sqrt(fabs(e)) * sqrt(fabs(c));
    } else {
        e = -(a / b * (c / b)) + 1.0f64;
        d = sqrt(fabs(e)) * fabs(b);
    }
    if e >= 0.0 {
        d = if b >= 0.0 { -d } else { d };
        *lr = (-b + d) / a;
        *sr = if *lr != 0.0 { c / *lr / a } else { *sr };
    } else {
        *sr = -(b / a);
        *lr = *sr;
        *si = fabs(d / a);
        *li = -*si;
    };
}

#[test]
fn test_rpoly() {
    const MAX_REAL_ROOT_NUM: usize = 5;
    const MAX_COMPLEX_PAIR_NUM: usize = 5;
    const MAX_ROOT_NUM: usize = MAX_REAL_ROOT_NUM + 2 * MAX_COMPLEX_PAIR_NUM;
    const MDP1: usize = MAX_ROOT_NUM + 1;

    use rand::{thread_rng, Rng};

    let mut rng = thread_rng();

    for _ in 0..1 {
        let real_root_num = rng.gen_range(0..=MAX_REAL_ROOT_NUM);
        let complex_root_pair_num = rng.gen_range(0..=MAX_COMPLEX_PAIR_NUM);

        let mut real_roots = [0.0; MAX_REAL_ROOT_NUM];
        for i in 0..real_root_num {
            if rng.gen_bool(0.9) {
                real_roots[i] = rng.gen_range(-10000.0..10000.0);
            }
        }

        let mut complex_roots_re = [0.0; MAX_COMPLEX_PAIR_NUM];
        let mut complex_roots_im = [0.0; MAX_COMPLEX_PAIR_NUM];
        for i in 0..complex_root_pair_num {
            complex_roots_im[i] = rng.gen_range(-10000.0..10000.0);
            if rng.gen_bool(0.9) {
                complex_roots_re[i] = rng.gen_range(-10000.0..10000.0);
            }
        }

        // generate coefficients
        let mut op = [0.0; MDP1];
        op[0] = 1.0;
        for i in 0..real_root_num {
            let n = i + 1;
            for j in 1..n + 1 {
                op[j] += op[j - 1] * (-real_roots[i]);
            }
        }
        for i in 0..complex_root_pair_num {
            let n = real_root_num + 2 * i + 1;
            let a = complex_roots_re[i];
            let b = complex_roots_im[i];
            for j in 1..n + 2 {
                if j - 1 < n {
                    op[j] += op[j - 1] * (-2.0 * a);
                }
                if j >= 2 {
                    op[j] += op[j - 2] * (a * a + b * b);
                }
            }
        }

        let mut Degree = real_root_num + 2 * complex_root_pair_num;
        let mut zeror = [0.0; MDP1];
        let mut zeroi = [0.0; MDP1];
        rpoly_ak1(&op, &mut Degree, &mut zeror, &mut zeroi);
        dbg!(Degree);
        dbg!(zeror);
        dbg!(zeroi);
    }

    // let mut Degree = 4;
    // let op = [
    //     323752928628.60541,
    //     9730768832.4398994,
    //     128621978.7391808,
    //     0.039863987093239961,
    //     3.0887751478756638e-12,
    // ];
    // let mut zeror = [0.0; 5];
    // let mut zeroi = [0.0; 5];
    // rpoly_ak1(&op, &mut Degree, &mut zeror, &mut zeroi);
    // dbg!(Degree);
    // dbg!(zeror);
    // dbg!(zeroi);
}
