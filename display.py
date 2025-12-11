import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np


def display_discs(discs, qx, qa, qy):
    fig, ax = plt.subplots(figsize=(10, 8))

    # Draw parallelogram boundary
    vertices = [
        (0, 0),
        (qx, 0),
        (qx + qa, qy),
        (qa, qy),
        (0, 0)  # close the shape
    ]
    xs, ys = zip(*vertices)
    ax.plot(xs, ys, 'k-', linewidth=2, label='Cell boundary')

    # Draw discs
    colors = plt.cm.Set3(np.linspace(0, 1, len(discs)))
    for idx, (center, radius) in enumerate(discs):
        circle = patches.Circle(center, radius, fill=True, alpha=0.3,
                                edgecolor='black', linewidth=1.5,
                                facecolor=colors[idx])
        ax.add_patch(circle)
        # Mark center
        ax.plot(center[0], center[1], 'ko', markersize=4)

    from main import genIntPointsInCell
    int_points = genIntPointsInCell(qx, qa, qy)
    px, py = zip(*int_points)
    ax.plot(px, py, 'r.', markersize=8, label='Integer points')

    # Configure axes
    ax.set_aspect('equal')
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('y', fontsize=12)

    ax.grid(True, alpha=0.3)
    ax.set_title(f'Parallelogram Cell (qx={qx}, qa={qa}, qy={qy})', fontsize=14)
    ax.legend(loc='upper right')

    # Set plot limits with some padding
    margin = max(qx, qa, qy) * 0.1
    ax.set_xlim(-margin, qx + qa + margin)
    ax.set_ylim(-margin, qy + margin)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # Example usage
    example_discs = [((0,0),2.5),((5,0),2.5),((2,5),2.5),((7,5),2.5),((4.5,3),0.5),((2.5,2),0.5)]
    qx = 5
    qa = 2
    qy = 5
    display_discs(example_discs, qx,qa,qy)
