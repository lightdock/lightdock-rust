// Default random number generator seed
pub const DEFAULT_SEED: u64 = 324_324;

// Translation interpolation step
pub const DEFAULT_TRANSLATION_STEP: f64 = 0.5;

// Rotation interpolation step
pub const DEFAULT_ROTATION_STEP: f64 = 0.5;

// When a quaternion SLERP is considered linear and not spherical
pub const LINEAR_THRESHOLD: f64 = 0.9995;

// Atomic contact is below this value
pub const INTERFACE_CUTOFF: f64 = 3.9;
pub const INTERFACE_CUTOFF2: f64 = INTERFACE_CUTOFF * INTERFACE_CUTOFF;

// Parsed PDB structures by lightdock start with this prefix
pub const DEFAULT_LIGHTDOCK_PREFIX: &str = "lightdock_";

// Membrane penalty for biasing the scoring
pub const MEMBRANE_PENALTY_SCORE: f64 = 999.0;

// ANM interpolation step
pub const DEFAULT_NMODES_STEP: f64 = 0.5;

// 1D NumPy arrays containing calculated ANM from ProDy
pub const DEFAULT_REC_NM_FILE: &str = "rec_nm.npy";
pub const DEFAULT_LIG_NM_FILE: &str = "lig_nm.npy";