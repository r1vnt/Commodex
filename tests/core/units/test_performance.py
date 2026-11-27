"""Opt-in correctness-checked timings, throughput and memory observations."""
from concurrent.futures import ThreadPoolExecutor
from decimal import Decimal as D
import json
import multiprocessing
import os
from pathlib import Path
import subprocess
import sys
import time
import tracemalloc
import pytest
import core.units as u
from core.units.batch import convert_many

pytestmark=pytest.mark.performance


def _worker(value): return u.convert(u.Q_(value,'kg'),'g').magnitude


def _pool():
    with multiprocessing.get_context('spawn').Pool(2) as pool: return pool.map(_worker,range(20))


def _threaded():
    with ThreadPoolExecutor(max_workers=2) as executor:
        return list(executor.map(_worker,range(20)))


def _operation(name):
    return {
        'warm_registry':lambda:u.get_registry(),
        'unit_parsing':lambda:u.parse_unit('USD/kg'),
        'quote_parsing':lambda:u.parse_quote('12.5 USD/kg'),
        'scalar_conversion':lambda:u.convert(u.Q_(D(10),'kg'),'g'),
        'batch_conversion':lambda:convert_many([u.Q_(D(1),'kg')]*100,'g'),
        'decimal_arithmetic':lambda:sum((D('.1') for _ in range(1000)),D(0)),
        'assay_construction':lambda:u.CrudeAssay.brent(),
        'context_creation':lambda:u.CrudeAssay.brent().as_context(),
        'contextual_conversion':lambda:u.convert(u.Q_(D(100),'bbl'),'kg',assay=u.CrudeAssay.brent()),
        'fx_triangulation':lambda:u.triangulate_fx_rate({'USD/EUR':D('.9'),'EUR/GBP':D('.8')},'USD','GBP'),
        'thermodynamic_solver':lambda:u.calculate_z_factor(u.Q_(1000,'kPa'),u.Q_(300,'K'),method='dak'),
        'serialisation':lambda:u.quantity_to_dict(u.Q_(D('123.45'),'USD/MWh_e'),include_metadata=True),
        'cache':lambda:u.decompose_quote('USD/kg'),
        'concurrent_conversion':_threaded,
        'multiprocessing':_pool,
    }[name]


@pytest.mark.parametrize('name',['warm_registry','unit_parsing','quote_parsing','scalar_conversion','batch_conversion','decimal_arithmetic','assay_construction','context_creation','contextual_conversion','fx_triangulation','thermodynamic_solver','serialisation','cache','concurrent_conversion','multiprocessing'])
def test_operation_benchmarks(name,request,record_property):
    pytest.importorskip('pytest_benchmark')
    benchmark=request.getfixturevalue('benchmark')
    result=benchmark.pedantic(_operation(name),rounds=3,iterations=1,warmup_rounds=0)
    assert result is not None
    assert benchmark.stats.stats.mean > 0
    record_property('operation',name)


def test_cold_registry_startup_benchmark(request):
    pytest.importorskip('pytest_benchmark')
    benchmark=request.getfixturevalue('benchmark')
    result=benchmark.pedantic(lambda:subprocess.run([sys.executable,'-B','-c','import core.units as u; u.verify_registry_integrity()'],capture_output=True,text=True,timeout=30),rounds=3,iterations=1)
    assert result.returncode == 0,result.stderr


def test_million_observation_throughput_and_dtype(record_property):
    np=pytest.importorskip('numpy')
    data=np.arange(1_000_000,dtype=np.float64)
    started=time.perf_counter()
    result=u.Q_(data,'kg').to('g')
    elapsed=time.perf_counter()-started
    np.testing.assert_allclose(result.magnitude,data*1000,rtol=1e-14)
    assert result.shape == (1_000_000,)
    record_property('observations_per_second',1_000_000/max(elapsed,1e-12))


def test_sustained_load_and_memory_consumption(record_property):
    tracemalloc.start()
    try:
        started=time.perf_counter()
        for _ in range(10000): assert _worker(1) == 1000
        current,peak=tracemalloc.get_traced_memory()
    finally: tracemalloc.stop()
    assert peak < 64*1024*1024
    record_property('elapsed_seconds',time.perf_counter()-started)
    record_property('peak_traced_bytes',peak)


def test_reviewed_performance_regression_baseline():
    path=os.environ.get('COMMODEX_BENCHMARK_BASELINE')
    if not path: pytest.skip('No reviewed same-host benchmark baseline supplied; timings are recorded, not a universal speed guarantee')
    baseline=json.loads(Path(path).read_text(encoding='utf-8'))
    assert baseline['python'] == sys.version.split()[0], 'Benchmark baseline interpreter mismatch'
    for name,limits in baseline['operations'].items():
        samples=[]
        for _ in range(5):
            start=time.perf_counter(); _operation(name)(); samples.append(time.perf_counter()-start)
        assert min(samples) <= limits['seconds']*limits.get('maximum_ratio',2)
