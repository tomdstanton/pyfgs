fn main() -> pyo3_stub_gen::Result<()> {
    let info = _pyfgs::stub_info()?;
    info.generate()?;
    Ok(())
}
