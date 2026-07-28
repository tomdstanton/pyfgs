pub mod models;
pub mod readers;
pub mod writers;
pub mod hmm;
pub mod cli;

use pyo3::prelude::*;
use pyo3_stub_gen::define_stub_info_gatherer;

// Re-exports
pub use models::*;
pub use readers::*;
pub use writers::*;
pub use hmm::*;

#[pymodule]
fn _pyfgs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Model>()?;
    m.add_class::<FastaReader>()?;
    m.add_class::<FastqReader>()?;
    m.add_class::<Gene>()?;
    m.add_class::<GeneBatch>()?;
    m.add_class::<GeneFinder>()?;
    m.add_class::<Mutation>()?;
    m.add_class::<BedWriter>()?;
    m.add_class::<VcfWriter>()?;
    m.add_class::<Gff3Writer>()?;
    m.add_class::<FnaWriter>()?;
    m.add_class::<FaaWriter>()?;
    
    // HMM parameter models
    m.add_class::<HmmGlobal>()?;
    m.add_class::<HmmLocal>()?;
    Ok(())
}

define_stub_info_gatherer!(stub_info);
