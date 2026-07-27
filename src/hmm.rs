use pyo3::prelude::*;
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};
use numpy::{IntoPyArray, ToPyArray, PyArray1, PyArray2, PyArray3};
use frag_gene_scan_rs::hmm::{Global, Local, PERIOD};

// Constants from frag_gene_scan_rs
const ACGT: usize = 4;
const BI_ACGT: usize = 16;
const TRI_ACGT: usize = 64;
const WINDOW: usize = 61;

/// A wrapper around the FragGeneScanRs Global HMM states
#[gen_stub_pyclass]
#[pyclass(unsendable)]
pub struct HmmGlobal {
    // We hold a pointer to the global state from GeneFinder
    global: *const Global,
}

impl HmmGlobal {
    pub fn new(global: &Global) -> Self {
        Self { global: global as *const Global }
    }
    
    fn get(&self) -> &Global {
        unsafe { &*self.global }
    }
}

#[gen_stub_pymethods]
#[pymethods]
impl HmmGlobal {
    /// Initial state probabilities
    #[pyo3(signature = ())]
    fn pi<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.get().pi.to_pyarray(py)
    }

    /// Transition probabilities (MM, MI, MD, II, IM, DD, DM, GE, GG, ER, RS, RR, ES, ES1)
    #[pyo3(signature = ())]
    fn transitions<'py>(&self, py: Python<'py>) -> pyo3::PyResult<Bound<'py, pyo3::types::PyDict>> {
        let dict = pyo3::types::PyDict::new(py);
        let tr = &self.get().tr;
        dict.set_item("MM", tr.mm)?;
        dict.set_item("MI", tr.mi)?;
        dict.set_item("MD", tr.md)?;
        dict.set_item("II", tr.ii)?;
        dict.set_item("IM", tr.im)?;
        dict.set_item("DD", tr.dd)?;
        dict.set_item("DM", tr.dm)?;
        dict.set_item("GE", tr.ge)?;
        dict.set_item("GG", tr.gg)?;
        dict.set_item("ER", tr.er)?;
        dict.set_item("RS", tr.rs)?;
        dict.set_item("RR", tr.rr)?;
        dict.set_item("ES", tr.es)?;
        dict.set_item("ES1", tr.es1)?;
        Ok(dict)
    }

    /// Insertion-to-insertion transition matrix [4 x 4]
    #[pyo3(signature = ())]
    fn tr_ii<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        let mut flat = Vec::with_capacity(ACGT * ACGT);
        for row in &self.get().tr_ii {
            flat.extend_from_slice(row);
        }
        numpy::ndarray::Array2::from_shape_vec([ACGT, ACGT], flat).unwrap().into_pyarray(py)
    }

    /// Match-to-insertion transition matrix [4 x 4]
    #[pyo3(signature = ())]
    fn tr_mi<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        let mut flat = Vec::with_capacity(ACGT * ACGT);
        for row in &self.get().tr_mi {
            flat.extend_from_slice(row);
        }
        numpy::ndarray::Array2::from_shape_vec([ACGT, ACGT], flat).unwrap().into_pyarray(py)
    }
}

/// A wrapper around the FragGeneScanRs Local (GC-specific) HMM states
#[gen_stub_pyclass]
#[pyclass(unsendable)]
pub struct HmmLocal {
    local: *const Local, // We use a raw pointer because locals are owned by GeneFinder
}

impl HmmLocal {
    pub fn new(local: &Local) -> Self {
        Self { local: local as *const Local }
    }
    
    fn get(&self) -> &Local {
        unsafe { &*self.local }
    }
}

#[gen_stub_pymethods]
#[pymethods]
impl HmmLocal {
    /// Emission probabilities for match states [PERIOD x BI_ACGT x ACGT]
    #[pyo3(signature = ())]
    fn e_m<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray3<f64>> {
        let local = self.get();
        let mut flat = Vec::with_capacity(PERIOD * BI_ACGT * ACGT);
        for i in 0..PERIOD {
            for j in 0..BI_ACGT {
                flat.extend_from_slice(&local.e_m[i][j]);
            }
        }
        numpy::ndarray::Array3::from_shape_vec([PERIOD, BI_ACGT, ACGT], flat).unwrap().into_pyarray(py)
    }

    /// Emission probabilities for match reverse states [PERIOD x BI_ACGT x ACGT]
    #[pyo3(signature = ())]
    fn e_m1<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray3<f64>> {
        let local = self.get();
        let mut flat = Vec::with_capacity(PERIOD * BI_ACGT * ACGT);
        for i in 0..PERIOD {
            for j in 0..BI_ACGT {
                flat.extend_from_slice(&local.e_m1[i][j]);
            }
        }
        numpy::ndarray::Array3::from_shape_vec([PERIOD, BI_ACGT, ACGT], flat).unwrap().into_pyarray(py)
    }

    /// Background noncoding transition matrix [4 x 4]
    #[pyo3(signature = ())]
    fn tr_rr<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        let local = self.get();
        let mut flat = Vec::with_capacity(ACGT * ACGT);
        for row in &local.tr_rr {
            flat.extend_from_slice(row);
        }
        numpy::ndarray::Array2::from_shape_vec([ACGT, ACGT], flat).unwrap().into_pyarray(py)
    }

    /// Start state transitions [WINDOW x TRI_ACGT]
    #[pyo3(signature = ())]
    fn tr_s<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        let local = self.get();
        let mut flat = Vec::with_capacity(WINDOW * TRI_ACGT);
        for row in &local.tr_s {
            flat.extend_from_slice(row);
        }
        numpy::ndarray::Array2::from_shape_vec([WINDOW, TRI_ACGT], flat).unwrap().into_pyarray(py)
    }

    /// End state transitions [WINDOW x TRI_ACGT]
    #[pyo3(signature = ())]
    fn tr_e<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        let local = self.get();
        let mut flat = Vec::with_capacity(WINDOW * TRI_ACGT);
        for row in &local.tr_e {
            flat.extend_from_slice(row);
        }
        numpy::ndarray::Array2::from_shape_vec([WINDOW, TRI_ACGT], flat).unwrap().into_pyarray(py)
    }
    
    /// Reverse start state transitions [WINDOW x TRI_ACGT]
    #[pyo3(signature = ())]
    fn tr_s1<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        let local = self.get();
        let mut flat = Vec::with_capacity(WINDOW * TRI_ACGT);
        for row in &local.tr_s1 {
            flat.extend_from_slice(row);
        }
        numpy::ndarray::Array2::from_shape_vec([WINDOW, TRI_ACGT], flat).unwrap().into_pyarray(py)
    }

    /// Reverse end state transitions [WINDOW x TRI_ACGT]
    #[pyo3(signature = ())]
    fn tr_e1<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        let local = self.get();
        let mut flat = Vec::with_capacity(WINDOW * TRI_ACGT);
        for row in &local.tr_e1 {
            flat.extend_from_slice(row);
        }
        numpy::ndarray::Array2::from_shape_vec([WINDOW, TRI_ACGT], flat).unwrap().into_pyarray(py)
    }
}
