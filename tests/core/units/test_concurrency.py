"""Task/thread/process isolation, recovery and bounded resource checks."""
import asyncio
from concurrent.futures import ThreadPoolExecutor
from decimal import Decimal as D
import gc
import multiprocessing
import subprocess
import sys
import tracemalloc
import pytest
import core.units as u


def test_async_fx_context_isolation_nested_context_cleanup_and_exception_recovery():
    async def worker(rate,fail=False):
        context = u.create_fx_context({"USD":D(1),"EUR":D(rate)})
        with u.ureg.context(context):
            await asyncio.sleep(.001)
            assert u.Q_(D(100),"USD").to("EUR").magnitude == D(100)*D(rate)
            if fail: raise RuntimeError("declared failure")
        return rate
    async def run():
        results = await asyncio.gather(worker('.8'),worker('.9'),worker('.7',True),return_exceptions=True)
        assert isinstance(results[-1],RuntimeError)
    asyncio.run(run())
    assert not u.Q_(1,"USD").is_compatible_with("EUR")


def test_concurrent_assay_conversion_does_not_cross_contaminate():
    assays = [u.CrudeAssay.brent(),u.CrudeAssay.wti(),u.CrudeAssay.dubai()]
    expected = [u.convert(u.Q_(100,"bbl"),"mt",assay=assay).magnitude for assay in assays]
    def convert(index):
        return index,u.convert(u.Q_(100,"bbl"),"mt",assay=assays[index]).magnitude
    with ThreadPoolExecutor(max_workers=8) as executor:
        results = list(executor.map(convert,[0,1,2]*80))
    assert all(number == expected[index] for index,number in results)


def test_registry_initialisation_and_lock_order_do_not_deadlock():
    program = """
from concurrent.futures import ThreadPoolExecutor
import core.units as u
with ThreadPoolExecutor(max_workers=8) as executor:
    values=list(executor.map(lambda _: (id(u.get_registry()),u.parse_unit('kg')),range(200)))
assert len({value[0] for value in values}) == 1
"""
    result = subprocess.run([sys.executable,"-B","-c",program],capture_output=True,text=True,timeout=25)
    assert result.returncode == 0,result.stderr


@pytest.mark.parametrize("method", ["spawn","fork","forkserver"])
def test_multiprocessing_start_method_preserves_registry_schema(method):
    if method not in multiprocessing.get_all_start_methods():
        pytest.skip(method+" is unavailable on this operating system; covered in the OS CI matrix")
    from tests.core.units.test_serialisation import _restore_in_child
    context = multiprocessing.get_context(method)
    parent,child = context.Pipe(duplex=False)
    payload = u.quantity_to_dict(u.Q_(D(2),"USD/MWh_e"),include_metadata=True)
    process = context.Process(target=_restore_in_child,args=(payload,child))
    process.start()
    try:
        assert parent.poll(20)
        assert parent.recv() == payload
    finally:
        process.join(20)
        if process.is_alive(): process.terminate(); process.join(5)
        parent.close(); child.close()
    assert process.exitcode == 0


def test_repeated_warm_conversions_do_not_retain_unbounded_temporary_quantities():
    for _ in range(200): u.convert(u.Q_(1,"kg"),"g")
    gc.collect()
    tracemalloc.start()
    try:
        before = tracemalloc.get_traced_memory()[0]
        for _ in range(2000): u.convert(u.Q_(1,"kg"),"g")
        gc.collect()
        retained = tracemalloc.get_traced_memory()[0]-before
    finally: tracemalloc.stop()
    assert retained < 2_000_000
