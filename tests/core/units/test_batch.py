"""Batch ordering, failure atomicity and optional array/table storage contracts."""
from decimal import Decimal as D
import json
import pytest
import core.units as u
from core.units.batch import convert_many


def test_empty_mixed_unit_batches_and_scalar_equivalence():
    assert convert_many([],"kg") == ()
    values = [u.Q_(D(1),"kg"),u.Q_(D(1000),"g"),u.Q_(D(1),"mt")]
    assert convert_many(values,"kg") == tuple(u.convert(value,"kg") for value in values)
    assert [v.magnitude for v in convert_many(values,"kg")] == [1,1,1000]
    assert all(isinstance(value.magnitude,D) for value in convert_many(values,"kg"))


def test_mixed_currency_batch_uses_explicit_fx_and_no_partial_result_escapes():
    output = convert_many([u.Q_(D(10),"USD"),u.Q_(D(10),"EUR")],"EUR",rates={"USD/EUR":D(".9")})
    assert [value.magnitude for value in output] == [9,10]
    with pytest.raises((ValueError,u.MissingFXRateError)):
        convert_many([u.Q_(1,"USD"),u.Q_(1,"JPY")],"EUR",rates={"USD/EUR":D(".9")})


def test_vectorised_numpy_broadcast_shape_dtype_and_empty_arrays():
    np = pytest.importorskip("numpy")
    data = np.arange(12,dtype=np.float64).reshape(3,4)
    output = u.Q_(data,"kg").to("g")
    np.testing.assert_allclose(output.magnitude,data*1000)
    assert output.shape == data.shape and output.magnitude.dtype == np.float64
    np.testing.assert_allclose((u.Q_(data,"kg")+u.Q_(np.array([1,2,3,4]),"kg")).magnitude,data+np.array([1,2,3,4]))
    assert u.Q_(np.empty((0,3)),"kg").to("g").shape == (0,3)
    # The validated scalar dispatcher deliberately does not accept an array.
    with pytest.raises((TypeError,ValueError,ArithmeticError)):
        u.convert(u.Q_(data,"kg"),"g")


def test_numpy_float32_and_missing_value_policy_is_explicit():
    np = pytest.importorskip("numpy")
    value = u.Q_(np.array([1,2,np.nan],dtype=np.float32),"kg").to("g")
    assert value.magnitude.dtype == np.float32
    assert np.isnan(value.magnitude[-1])
    with pytest.raises((ValueError,TypeError)):
        convert_many([u.Q_(1,"kg"),u.Q_(float('nan'),"kg")],"g")


def test_pandas_unit_aware_series_conversion():
    pd = pytest.importorskip("pandas")
    pint_pandas = pytest.importorskip("pint_pandas")
    pint_pandas.PintType.ureg = u.ureg
    series = pd.Series([1.,2.,3.],dtype="pint[kg]")
    output = series.pint.to("g")
    assert output.pint.magnitude.tolist() == [1000.,2000.,3000.]
    assert list(output.index) == list(series.index)


def test_pyarrow_schema_preserves_decimal_units_and_metadata():
    pa = pytest.importorskip("pyarrow")
    field = pa.field("amount",pa.decimal128(20,4),metadata={b"unit":b"MWh_e"})
    table = pa.Table.from_arrays([pa.array([D("1.2500"),None],type=field.type)],schema=pa.schema([field]))
    sink = pa.BufferOutputStream()
    with pa.ipc.new_stream(sink,table.schema) as writer: writer.write_table(table)
    restored = pa.ipc.open_stream(sink.getvalue()).read_all()
    assert restored.schema.field("amount").metadata[b"unit"] == b"MWh_e"
    assert restored.column(0).to_pylist() == [D("1.2500"),None]
