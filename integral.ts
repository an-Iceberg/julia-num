function range(from: number, to: number, step: number = 1): number[] {
  return [...Array(Math.floor((to - from) / step) + 1)]
    .map((_, i) => from + i * step);
}

/**
 * Calculates the integral of `f` between `a` and `b` using
 * [Simpson's ⅜ rule](https://en.wikipedia.org/wiki/Simpson%27s_rule#Composite_Simpson's_3/8_rule)
 * , so using polynomials of 3ʳᵈ degree.
 */
function int_2(
  a: number,
  b: number,
  f: (x: number) => number,
  h: number = 1e-2,
): number {
  let n = Math.ceil((b - a) / h);
  const x = (i: number): number => a + (i * h);

  // deno-fmt-ignore
  return (1/3)*h*range(1, n/2)
    .map(i => f(x(2*i - 2)) + 4*f(x(2*i - 1)) + f(x(2*i)))
    .reduce((prev, current) => prev + current); // Sum
}

/**
 * Calculates the integral of `f` between `a` and `b` using
 * [Boole's rule](https://en.wikipedia.org/wiki/Finite_difference_coefficient)
 * , so using polynomials of 4ʳᵈ degree.
 */
function int_3(
  a: number,
  b: number,
  f: (x: number) => number,
  h: number = 1e-2,
): number {
  let n = Math.ceil((b - a) / h);
  const x = (i: number): number => a + (i * h);

  // deno-fmt-ignore
  return (3/8)*h*range(1, n/3)
    .map(i => f(x(3*i - 3)) + 3*f(x(3*i - 2)) + 3*f(x(3*i - 1)) + f(x(3*i)))
    .reduce((prev, current) => prev + current); // Sum
}

/**
 * Calculates the integral of `f` between `a` and `b` using
 * [Boole's rule](https://en.wikipedia.org/wiki/Finite_difference_coefficient)
 * , so using polynomials of 4ʳᵈ degree.
 */
function int_4(
  a: number,
  b: number,
  f: (x: number) => number,
  h: number = 1e-2,
): number {
  let n = Math.ceil((b - a) / h);
  const x = (i: number): number => a + (i * h);

  // deno-fmt-ignore
  return (2/45)*h*range(1, n/4)
    .map((i) => 7*f(x(4*i - 4)) + 32*f(x(4*i - 3)) + 12*f(x(4*i - 2)) + 32*f(x(4*i - 1)) + 7*f(x(4*i)))
    .reduce((prev, current) => prev + current); // Sum
}

{
  const f = (x: number): number => Math.cos(x);
  const F = (x: number): number => Math.sin(x);
  let a = 0.0;
  let b = 5.0;

  // deno-fmt-ignore
  console.log(`Simpson's ⅓: ε = ${Math.abs(int_2(a, b, f) - (F(b) - F(a))).toExponential(2)}`);
  // deno-fmt-ignore
  console.log(`Simpson's ⅜: ε = ${Math.abs(int_3(a, b, f) - (F(b) - F(a))).toExponential(2)}`);
  // deno-fmt-ignore
  console.log(`Boole's    : ε = ${Math.abs(int_4(a, b, f) - (F(b) - F(a))).toExponential(2)}`);

  console.log();

  a = 30.;
  b = 40;

  // deno-fmt-ignore
  console.log(`Simpson's ⅓: ε = ${Math.abs(int_2(a, b, f) - (F(b) - F(a))).toExponential(2)}`);
  // deno-fmt-ignore
  console.log(`Simpson's ⅜: ε = ${Math.abs(int_3(a, b, f) - (F(b) - F(a))).toExponential(2)}`);
  // deno-fmt-ignore
  console.log(`Boole's    : ε = ${Math.abs(int_4(a, b, f) - (F(b) - F(a))).toExponential(2)}`);
}

console.log();

{
  const f = (x: number): number => 1 / x;
  const F = (x: number): number => Math.log(x);
  let a = 1.0;
  let b = 5.0;

  // deno-fmt-ignore
  console.log(`Simpson's ⅓: ε = ${Math.abs(int_2(a, b, f) - (F(b) - F(a))).toExponential(2)}`);
  // deno-fmt-ignore
  console.log(`Simpson's ⅜: ε = ${Math.abs(int_3(a, b, f) - (F(b) - F(a))).toExponential(2)}`);
  // deno-fmt-ignore
  console.log(`Boole's    : ε = ${Math.abs(int_4(a, b, f) - (F(b) - F(a))).toExponential(2)}`);

  console.log();

  a = 100.0;
  b = 150.0;

  // deno-fmt-ignore
  console.log(`Simpson's ⅓: ε = ${Math.abs(int_2(a, b, f) - (F(b) - F(a))).toExponential(2)}`);
  // deno-fmt-ignore
  console.log(`Simpson's ⅜: ε = ${Math.abs(int_3(a, b, f) - (F(b) - F(a))).toExponential(2)}`);
  // deno-fmt-ignore
  console.log(`Boole's    : ε = ${Math.abs(int_4(a, b, f) - (F(b) - F(a))).toExponential(2)}`);
}
