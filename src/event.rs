use std::collections::BTreeMap;

use serde::{Serialize, Deserialize};

/// Scattering event
#[derive(Debug, PartialEq, Default, Clone, Serialize, Deserialize)]
pub struct Event {
    pub number: i32,
    pub mpi: i32,
    pub scale: f64,
    pub alpha_qcd: f64,
    pub alpha_qed: f64,
    pub signal_process_id: i32,
    pub signal_process_vertex: i32,
    pub random_states: Vec<i32>,
    pub weights: Vec<f64>,
    pub weight_names: Vec<String>,
    pub vertices: Vec<Vertex>,
    pub xs: CrossSection,
    pub pdf_info: PdfInfo,
    pub energy_unit: String,
    pub length_unit: String,
}

#[derive(Debug, PartialEq, Default, Clone, Serialize, Deserialize)]
pub struct Vertex {
    pub barcode: i32,
    pub status: i32,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub t: f64,
    pub weights: Vec<f64>,
    pub particles_in: Vec<Particle>,
    pub particles_out: Vec<Particle>,
}

#[derive(Debug, PartialEq, Default, Clone, Serialize, Deserialize)]
pub struct Particle {
    pub id: i32,
    pub p: FourVector,
    pub m: f64,
    pub status: i32,
    pub theta: f64,
    pub phi: f64,
    pub flows: BTreeMap<i32, i32>,
    pub end_vtx: i32,
}

#[derive(Debug, PartialEq, PartialOrd, Default, Copy, Clone, Serialize, Deserialize)]
pub struct FourVector(pub [f64; 4]);

impl FourVector {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn txyz(t: f64, x: f64, y: f64, z: f64) -> Self {
        FourVector { 0: [t, x, y, z] }
    }
}

impl std::ops::Index<usize> for FourVector {
    type Output = f64;

    fn index(&self, idx: usize) -> &Self::Output {
        &self.0[idx]
    }
}

impl std::ops::IndexMut<usize> for FourVector {
    fn index_mut(&mut self, idx: usize) -> &mut Self::Output {
        &mut self.0[idx]
    }
}

#[derive(Debug, PartialEq, PartialOrd, Default, Copy, Clone, Serialize, Deserialize)]
pub struct CrossSection {
    pub cross_section: f64,
    pub cross_section_error: f64,
}

#[derive(Debug, PartialEq, PartialOrd, Default, Copy, Clone, Serialize, Deserialize)]
pub struct PdfInfo {
    pub parton_id: [i32; 2],
    pub x: [f64; 2],
    pub scale: f64,
    pub xf: [f64; 2],
    pub pdf_id: [i32; 2],
}
