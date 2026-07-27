use pyo3::prelude::*;
use pyo3::types::PyBytes;
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

use std::fs::File;
use seq_io::fasta::{Reader, Record as FastaRecord};
use seq_io::fastq::{self, Record as FastqRecord};

/// A memory-efficient FASTA parser yielding (header, sequence).
#[gen_stub_pyclass]
#[pyclass(unsendable)]
pub struct FastaReader {
    reader: Reader<File>,
}

#[gen_stub_pymethods]
#[pymethods]
impl FastaReader {
    /// Open a FASTA file for reading.
    ///
    /// Args:
    ///     path (str): The path to the FASTA file.
    #[new]
    fn new(path: String) -> PyResult<Self> {
        let file = File::open(&path)
            .map_err(|e| pyo3::exceptions::PyFileNotFoundError::new_err(format!("Could not open {}: {}", path, e)))?;
        Ok(Self { reader: Reader::new(file) })
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    // Pass `py: Python<'py>` to interact with Python memory
    fn __next__<'py>(mut slf: PyRefMut<'_, Self>, py: Python<'py>) -> PyResult<Option<(Bound<'py, PyBytes>, Bound<'py, PyBytes>)>> {
        match slf.reader.next() {
            Some(Ok(record)) => {
                let head = PyBytes::new(py, record.id_bytes());
                let seq = PyBytes::new(py, record.seq());
                Ok(Some((head, seq)))
            }
            Some(Err(e)) => Err(pyo3::exceptions::PyIOError::new_err(e.to_string())),
            None => Ok(None),
        }
    }
}

/// A memory-efficient FASTQ parser yielding (header, sequence, qualities).
#[gen_stub_pyclass]
#[pyclass(unsendable)]
pub struct FastqReader {
    reader: fastq::Reader<File>,
}

#[gen_stub_pymethods]
#[pymethods]
impl FastqReader {
    /// Open a FASTQ file for reading.
    #[new]
    fn new(path: String) -> PyResult<Self> {
        let file = File::open(&path)
            .map_err(|e| pyo3::exceptions::PyFileNotFoundError::new_err(format!("Could not open {}: {}", path, e)))?;
        Ok(Self { reader: fastq::Reader::new(file) })
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    /// Yields a tuple of (header, sequence, quality_string)
    fn __next__<'py>(mut slf: PyRefMut<'_, Self>, py: Python<'py>) -> PyResult<Option<(Bound<'py, PyBytes>, Bound<'py, PyBytes>, Bound<'py, PyBytes>)>> {
        match slf.reader.next() {
            Some(Ok(record)) => {
                let head = PyBytes::new(py, record.id_bytes());
                let seq = PyBytes::new(py, record.seq());
                let qual = PyBytes::new(py, record.qual());
                Ok(Some((head, seq, qual)))
            }
            Some(Err(e)) => Err(pyo3::exceptions::PyIOError::new_err(e.to_string())),
            None => Ok(None),
        }
    }
}
