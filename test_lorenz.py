#!/usr/bin/env python3
"""
Test script for Lorenz Attractor implementation
Tests the core functionality without GUI
"""

import numpy as np
from lorenz_attractor import LorenzAttractor


def test_lorenz_initialization():
    """Test Lorenz Attractor initialization"""
    print("Testing Lorenz Attractor initialization...")
    
    # Test default parameters
    lorenz = LorenzAttractor()
    assert lorenz.sigma == 10.0, "Default sigma should be 10.0"
    assert lorenz.rho == 28.0, "Default rho should be 28.0"
    assert abs(lorenz.beta - 8.0/3.0) < 1e-10, "Default beta should be 8/3"
    print("✓ Default initialization passed")
    
    # Test custom parameters
    lorenz = LorenzAttractor(sigma=5.0, rho=20.0, beta=2.0)
    assert lorenz.sigma == 5.0, "Custom sigma should be 5.0"
    assert lorenz.rho == 20.0, "Custom rho should be 20.0"
    assert lorenz.beta == 2.0, "Custom beta should be 2.0"
    print("✓ Custom initialization passed")


def test_lorenz_system():
    """Test Lorenz system equations"""
    print("\nTesting Lorenz system equations...")
    
    lorenz = LorenzAttractor()
    state = [1.0, 1.0, 1.0]
    t = 0.0
    
    derivatives = lorenz.lorenz_system(t, state)
    
    # Check that we get three derivatives
    assert len(derivatives) == 3, "Should return 3 derivatives"
    
    # Manually calculate expected values
    x, y, z = state
    expected_dx = lorenz.sigma * (y - x)  # 10 * (1 - 1) = 0
    expected_dy = x * (lorenz.rho - z) - y  # 1 * (28 - 1) - 1 = 26
    expected_dz = x * y - lorenz.beta * z  # 1 * 1 - (8/3) * 1 = 1 - 8/3 = -5/3
    
    assert abs(derivatives[0] - expected_dx) < 1e-10, f"dx/dt calculation error: {derivatives[0]} vs {expected_dx}"
    assert abs(derivatives[1] - expected_dy) < 1e-10, f"dy/dt calculation error: {derivatives[1]} vs {expected_dy}"
    assert abs(derivatives[2] - expected_dz) < 1e-10, f"dz/dt calculation error: {derivatives[2]} vs {expected_dz}"
    
    print("✓ Lorenz system equations passed")


def test_solve():
    """Test solving the Lorenz system"""
    print("\nTesting Lorenz system solver...")
    
    lorenz = LorenzAttractor()
    initial_state = [1.0, 1.0, 1.0]
    
    t, solution = lorenz.solve(initial_state, t_span=(0, 10), t_points=1000)
    
    # Check output shapes
    assert len(t) == 1000, "Should have 1000 time points"
    assert solution.shape == (1000, 3), "Solution should have shape (1000, 3)"
    
    # Check initial conditions
    assert abs(solution[0, 0] - initial_state[0]) < 1e-6, "Initial x should match"
    assert abs(solution[0, 1] - initial_state[1]) < 1e-6, "Initial y should match"
    assert abs(solution[0, 2] - initial_state[2]) < 1e-6, "Initial z should match"
    
    # Check that the solution evolves (values change)
    assert not np.allclose(solution[0], solution[-1]), "Solution should evolve over time"
    
    print("✓ Solver passed")
    print(f"  Initial state: {solution[0]}")
    print(f"  Final state: {solution[-1]}")


def test_different_initial_conditions():
    """Test with different initial conditions"""
    print("\nTesting different initial conditions...")
    
    lorenz = LorenzAttractor()
    
    # Test multiple initial conditions
    initial_conditions = [
        [0.0, 1.0, 0.0],
        [1.0, 1.0, 1.0],
        [5.0, 5.0, 5.0],
        [-1.0, -1.0, -1.0]
    ]
    
    for ic in initial_conditions:
        t, solution = lorenz.solve(ic, t_span=(0, 5), t_points=500)
        assert solution.shape == (500, 3), f"Solution shape incorrect for IC {ic}"
        assert abs(solution[0, 0] - ic[0]) < 1e-6, f"Initial x mismatch for IC {ic}"
        assert abs(solution[0, 1] - ic[1]) < 1e-6, f"Initial y mismatch for IC {ic}"
        assert abs(solution[0, 2] - ic[2]) < 1e-6, f"Initial z mismatch for IC {ic}"
        print(f"  ✓ IC {ic} passed")
    
    print("✓ Different initial conditions test passed")


def test_parameter_variation():
    """Test with different parameter values"""
    print("\nTesting parameter variations...")
    
    # Test different parameter sets
    parameter_sets = [
        {"sigma": 10.0, "rho": 28.0, "beta": 8.0/3.0},
        {"sigma": 5.0, "rho": 20.0, "beta": 2.0},
        {"sigma": 15.0, "rho": 35.0, "beta": 3.0}
    ]
    
    for params in parameter_sets:
        lorenz = LorenzAttractor(**params)
        t, solution = lorenz.solve([1.0, 1.0, 1.0], t_span=(0, 5), t_points=500)
        assert solution.shape == (500, 3), f"Solution shape incorrect for params {params}"
        print(f"  ✓ Parameters {params} passed")
    
    print("✓ Parameter variation test passed")


def main():
    """Run all tests"""
    print("=" * 60)
    print("Running Lorenz Attractor Tests")
    print("=" * 60)
    
    try:
        test_lorenz_initialization()
        test_lorenz_system()
        test_solve()
        test_different_initial_conditions()
        test_parameter_variation()
        
        print("\n" + "=" * 60)
        print("All tests passed! ✓")
        print("=" * 60)
        return 0
    except AssertionError as e:
        print(f"\n✗ Test failed: {e}")
        return 1
    except Exception as e:
        print(f"\n✗ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    exit(main())
