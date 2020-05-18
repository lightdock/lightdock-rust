use std::ops;
use std::f64;
use rand::Rng;
use std::f64::consts::PI;
use super::constants::LINEAR_THRESHOLD;

fn float_equals(x:f64, y:f64) -> bool {
    (x - y).abs() < f64::EPSILON
}


#[derive(Debug, Copy, Clone)]
pub struct Quaternion {
    pub w: f64,
    pub x: f64,
    pub y: f64,
    pub z: f64,
}


impl Quaternion {
    pub fn new(w:f64, x:f64, y:f64, z:f64) -> Quaternion {
        Quaternion {
            w,
            x,
            y,
            z
        }
    }

    pub fn conjugate(&self) -> Quaternion {
		Quaternion::new(self.w, -self.x, -self.y, -self.z)
	}

	pub fn dot(&self, other: Quaternion) -> f64 {
		self.w*other.w + self.x*other.x + self.y*other.y + self.z*other.z
	}

	pub fn norm2(&self) -> f64 {
		self.w*self.w + self.x*self.x + self.y*self.y + self.z*self.z
	}

	pub fn norm(&self) -> f64 {
		(self.w*self.w + self.x*self.x + self.y*self.y + self.z*self.z).sqrt()
	}

    pub fn normalize(&mut self) {
    	let norm = self.norm();
    	self.w /= norm;
        self.x /= norm;
        self.y /= norm;
        self.z /= norm;
    }

    pub fn inverse(&self) -> Quaternion{
    	self.conjugate() / self.norm2()
    }

    pub fn distance(&self, other: Quaternion) -> f64 {
    	let dot = self.dot(other);
    	1.0 - dot*dot
    }

    pub fn rotate(&self, vec3: Vec<f64>) -> Vec<f64> {
        let v = Quaternion::new(0., vec3[0], vec3[1], vec3[2]);
        let r = *self * v * self.inverse();
        vec![r.x, r.y, r.z]
    }

    pub fn lerp(&self, other: Quaternion, t: f64) -> Quaternion {
        *self * (1.0 - t) + other * t
    }

	pub fn slerp(&self, other: &Quaternion, t: f64) -> Quaternion {
		let mut q1 = *self;
		let mut q2 = *other;
		q1.normalize();
        q2.normalize();
        let mut q_dot = q1.dot(q2);

        // Patch to avoid the long path
        if q_dot < 0.0 {
            q1 = -q1;
            q_dot *= -1.0;
        }
        
        if q_dot > LINEAR_THRESHOLD {
            // Linear interpolation if quaternions are too close
            let mut result = q1 + (q2-q1) * t;
            result.normalize();
            result
        } else {
        	q_dot =((q_dot).min(1.0)).max(-1.0);
            let omega = q_dot.acos();
            let so = omega.sin();
            q1 * (((1.0-t)*omega).sin() / so) + q2 * ((t*omega).sin()/so)
        }
    }

    pub fn random(rng: &mut rand::prelude::StdRng) -> Quaternion {
        let u1 = rng.gen::<f64>();
        let u2 = rng.gen::<f64>();
        let u3 = rng.gen::<f64>();
        Quaternion::new(
            (1.0-u1).sqrt() * (2.0 * PI * u2).sin(),
            (1.0-u1).sqrt() * (2.0 * PI * u2).cos(),
            u1.sqrt() * (2.0 * PI * u3).sin(),
            u1.sqrt() * (2.0 * PI * u3).cos()
        )
    }
 }

 impl Default for Quaternion {
    fn default() -> Quaternion {
        Quaternion {
            w: 1.0,
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }
}

impl ops::Sub for Quaternion {
    type Output = Self;

    fn sub(self, other: Quaternion) -> Self::Output {
        Quaternion::new(self.w-other.w, self.x-other.x, self.y-other.y, self.z-other.z)
    }
}

impl ops::Add for Quaternion {
    type Output = Self;

    fn add(self, other: Quaternion) -> Self::Output {
        Quaternion::new(self.w+other.w, self.x+other.x, self.y+other.y, self.z+other.z)
    }
}

impl PartialEq for Quaternion {
    fn eq(&self, other: &Self) -> bool {
        float_equals(self.w, other.w) &&
        float_equals(self.x, other.x) &&
        float_equals(self.y, other.y) &&
        float_equals(self.z, other.z)
    }
}
impl Eq for Quaternion {}

impl ops::Neg for Quaternion {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Quaternion::new(-self.w, -self.x, -self.y, -self.z)
    }
}

impl ops::Mul<f64> for Quaternion {
    type Output = Self;

    fn mul(self, scalar: f64) -> Self::Output {
        Quaternion::new(scalar*self.w, scalar*self.x, scalar*self.y, scalar*self.z)
    }
}

impl ops::Mul for Quaternion {
	type Output = Self;

	fn mul(self, other: Quaternion) -> Self::Output {
		Quaternion::new(
			self.w * other.w - self.x * other.x - self.y * other.y - self.z * other.z,
        	self.w * other.x + self.x * other.w + self.y * other.z - self.z * other.y,
        	self.w * other.y - self.x * other.z + self.y * other.w + self.z * other.x,
        	self.w * other.z + self.x * other.y - self.y * other.x + self.z * other.w
        	)
	}
}

impl ops::Div<f64> for Quaternion {
    type Output = Self;

    fn div(self, scalar: f64) -> Self::Output {
        Quaternion::new(self.w/scalar, self.x/scalar, self.y/scalar, self.z/scalar)
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quaternion_default() {
        let q1:Quaternion = Default::default();
        assert_eq!(q1.w, 1.0);
        assert_eq!(q1.x, 0.0);
        assert_eq!(q1.y, 0.0);
        assert_eq!(q1.z, 0.0);
    }

    #[test]
    fn quaternion_sub() {
        let q1 = Quaternion::new(2.0, 0.0, 2.0, 0.0);
        let q2 = Quaternion::new(1.0, 0.0, 2.0, 1.0);
        let q3 = q1 - q2;
        assert_eq!(q3.w, 1.0);
        assert_eq!(q3.x, 0.0);
        assert_eq!(q3.y, 0.0);
        assert_eq!(q3.z, -1.0);
    }

    #[test]
    fn quaternion_add() {
        let q1 = Quaternion::new(2.0, -1.0, 2.0, 0.0);
        let q2 = Quaternion::new(1.0, 0.0, 2.0, 1.0);
        let q3 = q1 + q2;
        assert_eq!(q3.w, 3.0);
        assert_eq!(q3.x, -1.0);
        assert_eq!(q3.y, 4.0);
        assert_eq!(q3.z, 1.0);
    }

    #[test]
    fn quaternion_eq() {
    	let q1:Quaternion = Default::default();
    	let q2:Quaternion = Default::default();
        assert!(q1 == q2);
        let q3:Quaternion = Default::default();
    	let q4:Quaternion = Quaternion::new(1.000000000000001, 0.0, 0.0, 0.0);
        assert!(q3 != q4);
        let q5:Quaternion = Default::default();
    	let q6:Quaternion = Quaternion::new(1.0000000000000001, 0.0, 0.0, 0.0);
        assert!(q5 == q6);
    }

    #[test]
    fn quaternion_neg() {
        let q1 = Quaternion::new(2.0, -1.0, 2.0, 0.0);
        let q2 = Quaternion::new(-2.0, 1.0, -2.0, 0.0);
        assert!(q2 == -q1);
    }

    #[test]
    fn quaternion_mul_scalar() {
        let q1 = Quaternion::new(2.0, -1.0, 2.0, 0.0);
        let q2 = Quaternion::new(1.0, -0.5, 1.0, 0.0);
        assert!(q2 == q1 * 0.5);
    }

    #[test]
    fn quaternion_conjugate() {
        let q1 = Quaternion::new(2.0, -1.0, 2.0, 0.0);
        let q2 = Quaternion::new(2.0, 1.0, -2.0, 0.0);
        assert!(q2 == q1.conjugate());
    }

    #[test]
    fn quaternion_mul() {
        let q1 = Quaternion::new(1.0, 0.0, 0.0, 2.0);
        let q2 = Quaternion::new(3.0, -1.0, 4.0, 3.0);
        let expected = Quaternion::new(-3.0, -9.0, 2.0, 9.0);
        assert!(expected == q1*q2);

        let q1 = Quaternion::new(1.0, 0.0, 0.0, 2.0);
        let q2 = Quaternion::new(3.0, -1.0, 4.0, 3.0);
        let expected = Quaternion::new(-3.0, 7.0, 6.0, 9.0);
        assert!(expected == q2*q1);

        let q1 = Quaternion::new(1.0, 0.0, 0.0, 2.0);
        let q2 = Quaternion::new(3.0, -1.0, 4.0, 3.0);
        let q3 = Quaternion::new(1.0/2.0, -3.0, 2.0, 9.0);
        let expected = Quaternion::new(-147.0/2.0, 97.0/2.0, -93.0, 19.0/2.0);
        assert!(expected == q2*q1*q3);
    }

    #[test]
    fn test_conjugate_and_multiplication() {
        let q1 = Quaternion::new(1.0, 0.0, 0.0, 2.0);
        let q2 = Quaternion::new(3.0, -1.0, 4.0, 3.0);
        let expected = Quaternion::new(35.0, 0.0, 0.0, 0.0);
        assert!( (q1*q2).conjugate() == q2.conjugate() * q1.conjugate() );
        assert!(expected == q2.conjugate()*q2);
    }

    #[test]
    fn test_dot_product() {
        let q = Quaternion::new(2_f64.sqrt()/2.0, 0.0, 2_f64.sqrt()/2.0, 0.0);
        assert_eq!(1.0000000000000002, q.dot(q));
    }

    #[test]
    fn test_norm() {
        let q1 = Quaternion::new(1.0, -3.0, 4.0, 3.0);
        let q2 = Quaternion::new(3.0, -1.0, 4.0, 3.0);
        assert_eq!(5.916079783099616, q1.norm());
        assert_eq!((q1*q2).norm(), q1.norm()*q2.norm());
    }

    #[test]
    fn test_normalize() {
        let mut q1 = Quaternion::new(1.0, -3.0, 4.0, 3.0);
        let expected = Quaternion::new(0.1690308509457033, -0.50709255283711, 0.6761234037828132, 0.50709255283711);
        q1.normalize();
        assert!(expected == q1);
    }

    #[test]
    fn test_inverse() {
        let q1 = Quaternion::new(1.0, 0.0, 0.0, 2.0);
        let q2 = Quaternion::new(3.0, -1.0, 4.0, 3.0);
        let expected = Quaternion::new(-3.0/175.0, 9.0/175.0, -2.0/175.0, -9.0/175.0);
        assert!(expected == (q1*q2).inverse());
    }

    #[test]
    fn test_distance_is_zero() {
        let q = Quaternion::new(0.707106781, 0.0, 0.707106781, 0.0);
        assert_eq!(0.0000000010552720919321246, q.distance(q));
    }

    #[test]
    fn test_distance_is_one() {
        let q1 = Quaternion::new(0.707106781, 0.0, 0.707106781, 0.0);
        let q2 = Quaternion::new(0.707106781, 0.0, -0.707106781, 0.0);
        assert_eq!(1.0, q1.distance(q2));
    }

    #[test]
    fn test_distance_is_half() {
        let q1 = Quaternion::new(0.707106781, 0.0, 0.707106781, 0.0);
        let q2 = Quaternion::new(0.0, 0.0, 1.0, 0.0);
        assert_eq!(0.5000000002638181, q1.distance(q2));
    }

	#[test]
    fn test_distance_composite_rotation() {
        let q1 = Quaternion::new(1.0, 0.0, 0.0, 0.0);
        let q2 = Quaternion::new(0.5, 0.5, 0.5, 0.5);
        assert_eq!(0.75, q1.distance(q2));
    }

    #[test]
    fn test_rotation() {
        let q = Quaternion::new(0.707106781, 0.0, 0.707106781, 0.0);
        let v: Vec<f64> = vec![1.0, 0.0, 0.0];

        let v2 = q.rotate(v);

        assert_eq!(0.0, v2[0]);
        assert_eq!(0.0, v2[1]);
        assert_eq!(-1.0, v2[2]);
    }

    #[test]
    fn test_lerp_t_0() {
        let q1 = Quaternion::new(1.0, 0.0, 0.0, 2.0);
        let q2 = Quaternion::new(3.0, -1.0, 4.0, 3.0);

        let s = q1.lerp(q2, 0.0);

        assert!(s == q1);
    }

    #[test]
    fn test_lerp_t_1() {
        let q1 = Quaternion::new(1.0, 0.0, 0.0, 2.0);
        let q2 = Quaternion::new(3.0, -1.0, 4.0, 3.0);

        let s = q1.lerp(q2, 1.0);

        assert!(s == q2);
    }

    #[test]
    fn test_slerp_t_0() {
        let q1 = Quaternion::new(1.0, 0.0, 0.0, 2.0);
        let q2 = Quaternion::new(3.0, -1.0, 4.0, 3.0);
        let expected = Quaternion::new(0.4472135954999579, 0.0, 0.0, 0.8944271909999159);

        let s = q1.slerp(&q2, 0.0);

        assert!(expected == s);
    }

    #[test]
    fn test_slerp_t_1() {
        let q1 = Quaternion::new(1.0, 0.0, 0.0, 2.0);
        let q2 = Quaternion::new(3.0, -1.0, 4.0, 3.0);
        let expected = Quaternion::new(0.50709255283711, -0.1690308509457033, 0.6761234037828132, 0.50709255283711);

        let s = q1.slerp(&q2, 1.0);

        assert!(expected == s);
    }

    #[test]
    fn test_slerp_same_quaternion() {
        let q1 = Quaternion::new(0.7071067811865476, 0.0, 0.0, 0.7071067811865476);
        let q2 = Quaternion::new(0.7071067811865476, 0.0, 0.0, 0.7071067811865476);

        let s = q1.slerp(&q2, 0.1);

        assert!(s == q2);
    }

    #[test]
    fn test_slerp_t_half_y() {
        let q1 = Quaternion::new(1.0, 0.0, 0.0, 0.0);
        let q2 = Quaternion::new(0.0, 0.0, 1.0, 0.0);

        let s = q1.slerp(&q2, 0.5);

        let expected = Quaternion::new(0.7071067811865475, 0.0, 0.7071067811865475, 0.0);
        assert!(expected == s);
    }

    #[test]
    fn test_slerp_t_half_bank_zero() {
        let q1 = Quaternion::new(0.7071067811865475, 0.0, 0.0, 0.7071067811865475);
        let q2 = Quaternion::new(0.0, 0.7071067811865475, 0.7071067811865475, 0.0);

        let s = q1.slerp(&q2, 0.5);

        let expected = Quaternion::new(0.5, 0.5, 0.5, 0.5);
        assert!(expected == s);
    }

    #[test]
    fn test_random_quaternion() {
    	use rand::SeedableRng;
    	let mut rng = SeedableRng::seed_from_u64(324324324);
        let q = Quaternion::random(&mut rng);

        let expected = Quaternion::new(0.31924330894562036, -0.5980633213833059, 
        							   0.5444724265858514, 0.49391674399349367);
        assert!(expected == q);
    }
}